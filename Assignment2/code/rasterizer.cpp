// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f>& positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return { id };
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i>& indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return { id };
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f>& cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return { id };
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


/// <summary>
/// Implement this function to check if the point (x, y) is inside the triangle represented by v[0], v[1], v[2]
/// 判断点是否在三角形内
/// </summary>
static bool insideTriangle(float x, float y, const Vector3f* v, float winding)
{
    const int n = 3;
    // n 是顶点数量; i 表示当前顶点; j 表示上一个顶点
    // 灵魂代码是 j = i++
    for (int i = 0, j = n - 1; i < 3; j = i++) {
        // 顺时针时反转边方向，统一为内部位于边左侧的逆时针约定。
        float edgeX = winding * (v[i].x() - v[j].x());
        float edgeY = winding * (v[i].y() - v[j].y());
        float pX = x - v[j].x();
        float pY = y - v[j].y();
        float crossResult = edgeX * pY - edgeY * pX;
        // Y 向上、顶点逆时针：向下的边是左边，水平向左的边是上边。
        // 左上规则让共享边上的采样点只属于一个三角形。
        bool isTopLeft = edgeY < 0 || (edgeY == 0 && edgeX < 0);
        if (crossResult < 0 || (crossResult == 0 && !isTopLeft)) return false;
    }

    /* 简短形式
    for (int i = 0, j = n - 1; i < 3; j = i++) {
        Eigen::Vector2f edge(v[i].x() - v[j].x(), v[i].y() - v[j].y());
        Eigen::Vector2f ray(x - _v[j].x(), y - _v[j].y());
        auto res = edge.cross(ray);
        if (res < 0) return false;
    }
    */

    return true;
}

/// <summary>
/// 重心坐标
/// </summary>
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float doubleArea = (v[1].x() - v[0].x()) * (v[2].y() - v[0].y()) - (v[1].y() - v[0].y()) * (v[2].x() - v[0].x());
    float doubleAreaA = (v[1].x() - x) * (v[2].y() - y) - (v[1].y() - y) * (v[2].x() - x);
    float doubleAreaB = (v[2].x() - x) * (v[0].y() - y) - (v[2].y() - y) * (v[0].x() - x);
    float alpha = doubleAreaA / doubleArea;
    float beta = doubleAreaB / doubleArea;
    float gamma = 1 - alpha - beta;
    return { alpha, beta, gamma };
}

// 在选定的着色位置，使用原始裁剪空间 W 对顶点颜色进行透视校正插值。
static Eigen::Vector3f interpolateColor(float x, float y, const Triangle& t)
{
    auto [alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    float a = alpha / t.clip_w[0];
    float b = beta / t.clip_w[1];
    float c = gamma / t.clip_w[2];
    return (a * t.getColor(0) + b * t.getColor(1) + c * t.getColor(2)) / (a + b + c);
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    // 从 NDC [-1, 1] 映射到 [0.1, 50] （其实这一步在当前脚本没什么必要）
    float f1 = (50 - 0.1) / 2.0; // 缩放系数
    float f2 = (50 + 0.1) / 2.0; // 平移量

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };

        // Homogeneous division
        for (int k = 0; k < 3; ++k) {
            t.clip_w[k] = v[k].w();
            v[k] /= t.clip_w[k];
        }

        // Viewport transformation
        for (auto& vert : v)
        {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
        }

        auto color0 = col[i[0]];
        auto color1 = col[i[1]];
        auto color2 = col[i[2]];
        
        t.setColor(0, color0[0], color0[1], color0[2]);
        t.setColor(1, color1[0], color1[1], color1[2]);
        t.setColor(2, color2[0], color2[1], color2[2]);

        rasterize_triangle(t);
    }

    if (msaaEnable)
        resolve();
}

/// <summary>
/// Screen space rasterization
/// 屏幕空间光栅化
/// </summary>
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    // Cull Off：接收两种绕序，仅跳过没有面积的退化三角形。
    Eigen::Vector3f e1 = t.v[1] - t.v[0];
    Eigen::Vector3f e2 = t.v[2] - t.v[0];
    float area = e1.x() * e2.y() - e1.y() * e2.x();
    if (area == 0.0f) return;
    const float winding = area > 0.0f ? 1.0f : -1.0f;

    auto v = t.toVector4();

    int left = width, bottom = height;
    int right = 0, top = 0;
    for (auto& p : v) {
        left = MAX(MIN(left, floor(p.x())), 0);
        bottom = MAX(MIN(bottom, floor(p.y())), 0);
        right = MIN(MAX(right, ceil(p.x())), width);
        top = MIN(MAX(top, ceil(p.y())), height);
    }

    // 包围盒保持像素坐标，每个像素内部检查 2×2 个样本。
    const int samplesPerAxis = msaaEnable ? 2 : 1;
    for (int y = bottom; y < top; ++y) {
        for (int x = left; x < right; ++x) {
            float shadeX = x + 0.5f;
            float shadeY = y + 0.5f;
            bool shadingPositionChosen = false;
            bool shaded = false;
            Eigen::Vector3f pixelColor = Eigen::Vector3f::Zero();

            for (int j = 0; j < samplesPerAxis; ++j) {
                for (int i = 0; i < samplesPerAxis; ++i) {
                    float sampleX = x + (i + 0.5f) / samplesPerAxis;
                    float sampleY = y + (j + 0.5f) / samplesPerAxis;
                    if (!insideTriangle(sampleX, sampleY, t.v, winding)) continue;

                    // 质心采样：中心在三角形内就用中心，否则用第一个被覆盖的样本。
                    // 着色位置只依赖几何覆盖，在深度测试前确定。
                    if (!shadingPositionChosen) {
                        if (!insideTriangle(shadeX, shadeY, t.v, winding)) {
                            shadeX = sampleX;
                            shadeY = sampleY;
                        }
                        shadingPositionChosen = true;
                    }

                    // 深度仍在各个样本的位置独立插值和测试。
                    auto [alpha, beta, gamma] = computeBarycentric2D(sampleX, sampleY, t.v);
                    float depth = alpha * v[0].z() + beta * v[1].z() + gamma * v[2].z();
                    int sx = x * samplesPerAxis + i;
                    int sy = y * samplesPerAxis + j;
                    int index = get_index(sx, sy);
                    if (depth_buf[index] > depth) continue;

                    // 同一个三角形在同一个像素中最多着色一次，颜色由通过测试的样本共享。
                    if (!shaded) {
                        pixelColor = interpolateColor(shadeX, shadeY, t);
                        shaded = true;
                    }
                    depth_buf[index] = depth;
                    if (msaaEnable)
                        set_sample_pixel(Eigen::Vector3f(sx, sy, 0), pixelColor);
                    else
                        set_pixel(Eigen::Vector3f(x, y, 0), pixelColor);
                }
            }
        }
    }

    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
}

void rst::rasterizer::resolve()
{
    for (int y = 0; y < logicHeight; y += 2) {
        for (int x = 0; x < logicWidth; x += 2) {
            int index0 = get_index(x, y);
            int index1 = index0 + 1;
            int index2 = get_index(x, y + 1);
            int index3 = index2 + 1;
            auto targetColor = sample_color[index0] 
                + sample_color[index1]
                + sample_color[index2]
                + sample_color[index3];
            set_pixel(Eigen::Vector3f(x >> 1, y >> 1, 0), targetColor / 4);
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{ 0, 0, 0 });
        std::fill(sample_color.begin(), sample_color.end(), Eigen::Vector3f{ 0, 0, 0 });
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), -std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h, bool msaaEnable) : width(w), height(h)
{
    int baseSize = w * h;
    frame_buf.resize(baseSize);
    this->msaaEnable = msaaEnable;
    if (msaaEnable)
    {
        msaaScalar = 2;
        depth_buf.resize(baseSize * 4);
        sample_color.resize(baseSize * 4);
    }
    else
    {
        msaaScalar = 1;
        depth_buf.resize(baseSize);
        sample_color.reserve(0);
    }
    logicWidth = width * msaaScalar;
    logicHeight = height * msaaScalar;
}

int rst::rasterizer::get_index(int x, int y)
{
    return (logicHeight - 1 - y) * logicWidth + x;
}

/// <summary>
/// 光栅化时：左下角是原点，Y 向上增加。
/// 显示图像时：数组第一行是图像最上面的一行。
/// </summary>
/// <param name="point"> 表示像素位置，这里只用了 x、y，没有用 z </param>
/// <param name="color"> 表示 RGB 颜色 </param>
void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;
}

void rst::rasterizer::set_sample_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    auto ind = (logicHeight - 1 - point.y()) * logicWidth + point.x();
    sample_color[ind] = color;
}

// clang-format on
