//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_normals(const std::vector<Eigen::Vector3f>& normals)
{
    auto id = get_next_id();
    nor_buf.emplace(id, normals);

    normal_id = id;

    return {id};
}


// Bresenham's line drawing algorithm
void rst::rasterizer::draw_line(Eigen::Vector3f begin, Eigen::Vector3f end)
{
    auto x1 = begin.x();
    auto y1 = begin.y();
    auto x2 = end.x();
    auto y2 = end.y();

    Eigen::Vector3f line_color = {255, 255, 255};

    int x,y,dx,dy,dx1,dy1,px,py,xe,ye,i;

    dx=x2-x1;
    dy=y2-y1;
    dx1=fabs(dx);
    dy1=fabs(dy);
    px=2*dy1-dx1;
    py=2*dx1-dy1;

    if(dy1<=dx1)
    {
        if(dx>=0)
        {
            x=x1;
            y=y1;
            xe=x2;
        }
        else
        {
            x=x2;
            y=y2;
            xe=x1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point,line_color);
        for(i=0;x<xe;i++)
        {
            x=x+1;
            if(px<0)
            {
                px=px+2*dy1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    y=y+1;
                }
                else
                {
                    y=y-1;
                }
                px=px+2*(dy1-dx1);
            }
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point,line_color);
        }
    }
    else
    {
        if(dy>=0)
        {
            x=x1;
            y=y1;
            ye=y2;
        }
        else
        {
            x=x2;
            y=y2;
            ye=y1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point,line_color);
        for(i=0;y<ye;i++)
        {
            y=y+1;
            if(py<=0)
            {
                py=py+2*dx1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    x=x+1;
                }
                else
                {
                    x=x-1;
                }
                py=py+2*(dx1-dy1);
            }
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point,line_color);
        }
    }
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

static bool insideTriangle(float x, float y, const Vector4f* _v){
    Vector3f v[3];
    for(int i=0;i<3;i++)
        v[i] = {_v[i].x(),_v[i].y(), 1.0};
    // 两个点叉积，得到边的直线方程系数
    Vector3f f0,f1,f2;
    f0 = v[1].cross(v[0]);
    f1 = v[2].cross(v[1]);
    f2 = v[0].cross(v[2]);
    Vector3f p(x,y,1.);
    // 点积就是把点代入直线方程：判断目标点和三角形另一顶点是否在边的同一侧
    if((p.dot(f0)*f0.dot(v[2])>0) && (p.dot(f1)*f1.dot(v[0])>0) && (p.dot(f2)*f2.dot(v[1])>0))
        return true;
    return false;
}

/// <summary>
/// 重心坐标
/// </summary>
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector4f* v)
{
    float doubleArea = (v[1].x() - v[0].x()) * (v[2].y() - v[0].y()) - (v[1].y() - v[0].y()) * (v[2].x() - v[0].x());
    float doubleAreaA = (v[1].x() - x) * (v[2].y() - y) - (v[1].y() - y) * (v[2].x() - x);
    float doubleAreaB = (v[2].x() - x) * (v[0].y() - y) - (v[2].y() - y) * (v[0].x() - x);
    float alpha = doubleAreaA / doubleArea;
    float beta = doubleAreaB / doubleArea;
    float gamma = 1 - alpha - beta;
    return { alpha, beta, gamma };
}

/*
* 作业模板里的实现：容易引入误差，且有重复计算
* 类似这样的代码，有如下问题：v[1].x() * v[2].y() - v[2].x() * v[1].y()
* 屏幕坐标有几百，乘积就可能有几十万。但三角形很小时，最后求出的面积却很小。两个很接近的大数相减，容易损失有效数字，这叫消减误差（cancellation）。
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector4f* v){
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y())
        / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y())
        / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y())
        / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}
*/

void rst::rasterizer::draw(std::vector<Triangle *> &TriangleList) {

    // 从 NDC [-1, 1] 映射到 [0.1, 50] >> 近平面 1 映射到 50 while 远平面 -1 映射到 0.1
    // 浮点数在接近 0 时，能区分的数值间隔更小；把远处密集的深度放到接近 0 的位置，确实可能提高精度。 这也是反向深度（Reversed-Z）的核心思路
    // 反向深度的好处来自：把密集的部分搬到浮点数刻度更细的位置
    // 浮点表示方面： 远处密集的深度落到接近 0、浮点刻度更细的位置，更容易区分。
    // 插值运算方面： 使用这些较小的深度值计算，重心坐标误差带来的深度偏差也可能明显减小。
    // https://developer.nvidia.com/blog/visualizing-depth-precision/
    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (const auto& t:TriangleList)
    {
        Triangle newtri = *t;

        std::array<Eigen::Vector4f, 3> mm {
                (view * model * t->v[0]),
                (view * model * t->v[1]),
                (view * model * t->v[2])
        };

        std::array<Eigen::Vector3f, 3> viewspace_pos;

        std::transform(mm.begin(), mm.end(), viewspace_pos.begin(), [](auto& v) {
            return v.template head<3>();
        });

        Eigen::Vector4f v[] = {
                mvp * t->v[0],
                mvp * t->v[1],
                mvp * t->v[2]
        };

        // Homogeneous division
        // 与作业2不同，这里没有直接整个向量除以w，这样可以保留w，并且下面setVertex设置的也是4维的齐次坐标
        for (auto& vec : v) {
            vec.x()/=vec.w();
            vec.y()/=vec.w();
            vec.z()/=vec.w();
        }

        /*
        * 仿射变换保持直线和平行关系，但不保证保持角度
        * 所以如果我们要在观测空间计算光照，需要得到在观测空间真正垂直表面的法向量
        */
        Eigen::Matrix4f inv_trans = (view * model).inverse().transpose();
        Eigen::Vector4f n[] = {
                inv_trans * to_vec4(t->normal[0], 0.0f),
                inv_trans * to_vec4(t->normal[1], 0.0f),
                inv_trans * to_vec4(t->normal[2], 0.0f)
        };

        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            //screen space coordinates
            newtri.setVertex(i, v[i]);
        }

        for (int i = 0; i < 3; ++i)
        {
            //view space normal
            newtri.setNormal(i, n[i].head<3>());
        }

        newtri.setColor(0, 148,121.0,92.0);
        newtri.setColor(1, 148,121.0,92.0);
        newtri.setColor(2, 148,121.0,92.0);

        // Also pass view space vertice position
        rasterize_triangle(newtri, viewspace_pos);
    }
}

static Eigen::Vector3f interpolate(
    float alpha, float beta, float gamma,
    const Eigen::Vector3f& vert1,
    const Eigen::Vector3f& vert2,
    const Eigen::Vector3f& vert3,
    float weight)
{
    return (alpha * vert1 + beta * vert2 + gamma * vert3) / weight;
}

static Eigen::Vector2f interpolate(
    float alpha, float beta, float gamma,
    const Eigen::Vector2f& vert1,
    const Eigen::Vector2f& vert2,
    const Eigen::Vector2f& vert3,
    float weight)
{
    auto u = (alpha * vert1[0] + beta * vert2[0] + gamma * vert3[0]);
    auto v = (alpha * vert1[1] + beta * vert2[1] + gamma * vert3[1]);

    u /= weight;
    v /= weight;

    return Eigen::Vector2f(u, v);
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t, const std::array<Eigen::Vector3f, 3>& view_pos) 
{
    int left = width, right = 0, top = 0, bottom = height;
    for (const auto& p : t.v) {
        left = MAX(MIN(floor(p.x()), left), 0);
        right = MIN(MAX(ceil(p.x()), right), width);
        bottom = MAX(MIN(floor(p.y()), bottom), 0);
        top = MIN(MAX(ceil(p.y()), top), height);
    }

    for (int y = bottom; y < top; ++y) {
        for (int x = left; x < right; ++x) {
            float sampleX = x + 0.5, sampleY = y + 0.5;
            if (!insideTriangle(sampleX, sampleY, t.v))
                continue;

            auto [alpha, beta, gamma] = computeBarycentric2D(sampleX, sampleY, t.v);

            // 透视除法后的深度插值（仅用于深度测试是OK的）
            float depth = alpha * t.a().z() + beta * t.b().z() + gamma * t.c().z();
            int index = get_index(x, y);
            if (depth_buf[index] > depth)
                continue;
            depth_buf[index] = depth;

            float w1 = t.a().w(), w2 = t.b().w(), w3 = t.c().w();
            float w = 1 / (alpha / w1 + beta / w2 + gamma / w3);

            // 位置插值
            Vector3f position = view_pos[0] * alpha / w1 + view_pos[1] * beta / w2 + view_pos[2] * gamma / w3;
            position *= w;
            // 法向量插值
            Vector3f normal = t.normal[0] * alpha / w1 + t.normal[1] * beta / w2 + t.normal[2] * gamma / w3;
            normal *= w;
            normal.normalize();
            // 颜色插值
            Vector3f color = t.color[0] * alpha / w1 + t.color[1] * beta / w2 + t.color[2] * gamma / w3;
            color *= w;
            // 纹理插值
            Vector2f tex = t.tex_coords[0] * alpha / w1 + t.tex_coords[1] * beta / w2 + t.tex_coords[2] * gamma / w3;
            tex *= w;

            fragment_shader_payload fragment(color, normal, tex, t.tex);
            fragment.view_pos = position;
            color = fragment_shader(fragment);

            set_pixel(Eigen::Vector2i{ x, y }, color);
        }
    }


    // TODO: From your HW3, get the triangle rasterization code.
    // TODO: Inside your rasterization loop:
    //    * v[i].w() is the vertex view space depth value z.
    //    * Z is interpolated view space depth for the current pixel
    //    * zp is depth between zNear and zFar, used for z-buffer

    // float Z = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    // float zp = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    // zp *= Z;

    // TODO: Interpolate the attributes:
    // auto interpolated_color
    // auto interpolated_normal
    // auto interpolated_texcoords
    // auto interpolated_shadingcoords

    // Use: fragment_shader_payload payload( interpolated_color, interpolated_normal.normalized(), interpolated_texcoords, texture ? &*texture : nullptr);
    // Use: payload.view_pos = interpolated_shadingcoords;
    // Use: Instead of passing the triangle's color directly to the frame buffer, pass the color to the shaders first to get the final color;
    // Use: auto pixel_color = fragment_shader(payload);

 
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
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), -std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);

    texture = std::nullopt;
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-y)*width + x;
}

void rst::rasterizer::set_pixel(const Vector2i &point, const Eigen::Vector3f &color)
{
    //old index: auto ind = point.y() + point.x() * width;
    int ind = (height-point.y())*width + point.x();
    frame_buf[ind] = color;
}

void rst::rasterizer::set_vertex_shader(std::function<Eigen::Vector3f(vertex_shader_payload)> vert_shader)
{
    vertex_shader = vert_shader;
}

void rst::rasterizer::set_fragment_shader(std::function<Eigen::Vector3f(fragment_shader_payload)> frag_shader)
{
    fragment_shader = frag_shader;
}

