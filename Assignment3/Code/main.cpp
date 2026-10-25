#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <opencv2/opencv.hpp>

#include "global.hpp"
#include "rasterizer.hpp"
#include "Triangle.hpp"
#include "Shader.hpp"
#include "Texture.hpp"
#include "OBJ_Loader.h"

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float angle)
{
    Eigen::Matrix4f rotation;
    angle = angle * MY_PI / 180.f;
    rotation << cos(angle), 0, sin(angle), 0,
                0, 1, 0, 0,
                -sin(angle), 0, cos(angle), 0,
                0, 0, 0, 1;

    Eigen::Matrix4f scale;
    scale << 2.5, 0, 0, 0,
              0, 2.5, 0, 0,
              0, 0, 2.5, 0,
              0, 0, 0, 1;

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    return translate * rotation * scale;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    float t = tan(eye_fov / 2 / 180 * MY_PI) * abs(zNear);
    float b = -t;
    float r = t * aspect_ratio;
    float l = -r;

    // 课件里近平面和远平面都应该是负数
    float n = -abs(zNear);
    float f = -abs(zFar);

    // 透视到正交
    Eigen::Matrix4f perspToOrtho;
    perspToOrtho <<
        n, 0, 0, 0,
        0, n, 0, 0,
        0, 0, n + f, -n * f,
        0, 0, 1, 0;

    // 平移到原点
    Eigen::Matrix4f translate;
    translate <<
        1, 0, 0, -(l + r) / 2,
        0, 1, 0, -(t + b) / 2,
        0, 0, 1, -(n + f) / 2,
        0, 0, 0, 1;

    // 缩放到[-1, 1]
    Eigen::Matrix4f scale;
    scale <<
        2 / (r - l), 0, 0, 0,
        0, 2 / (t - b), 0, 0,
        0, 0, 2 / (n - f), 0,
        0, 0, 0, 1;

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    projection = scale * translate * perspToOrtho * projection;
    return projection;
}

Eigen::Vector3f vertex_shader(const vertex_shader_payload& payload)
{
    return payload.position;
}

Eigen::Vector3f normal_fragment_shader(const fragment_shader_payload& payload)
{
    // 法线可视化：把 [-1,1] 的法线映射到 [0,1] 后当作颜色
    Eigen::Vector3f return_color = (payload.normal.head<3>() + Eigen::Vector3f(1.0f, 1.0f, 1.0f)) / 2.f;
    return return_color * 255.f;
}

// ======================== 场景 / 材质 / 光照 ========================

struct light
{
    Eigen::Vector3f position;
    Eigen::Vector3f intensity;
};

// 一个材质需要的四项参数
// ka: 环境光反射系数; kd: 漫反射颜色; ks: 高光系数; p: 高光指数（越大高光越集中）
struct material
{
    Eigen::Vector3f ka;
    Eigen::Vector3f kd;
    Eigen::Vector3f ks;
    float p;
};

// 整个场景共享的参数：两盏点光源 + 环境光强度 + 相机位置
// 相机固定在 view space 的 (0,0,10)，和 main 里的 eye_pos 保持一致
struct scene
{
    std::vector<light> lights;
    Eigen::Vector3f amb;
    Eigen::Vector3f eye;
};

static const scene g_scene = {
    {
        light{ Eigen::Vector3f(20.f, 20.f, 20.f),   Eigen::Vector3f(500.f, 500.f, 500.f) },
        light{ Eigen::Vector3f(-20.f, 20.f, 0.f),   Eigen::Vector3f(500.f, 500.f, 500.f) }
    },
    Eigen::Vector3f(10.f, 10.f, 10.f),
    Eigen::Vector3f(0.f, 0.f, 10.f)
};

// ======================== 着色公共工具 ========================

// 纹理采样，返回 0~255 的颜色（和 Texture::getColor 的约定一致，不在这里归一化）
static Eigen::Vector3f sample_color(Texture* tex, const Eigen::Vector2f& uv, bool bilinear)
{
    if (tex == nullptr) return Eigen::Vector3f::Zero();
    return bilinear ? tex->getColorBilinear(uv.x(), uv.y())
                    : tex->getColor(uv.x(), uv.y());
}

// 把纹理颜色当成“高度”：取模长（灰度图三通道相同，等价于取通道值），归一到 [0,1]
static float sample_height(Texture* tex, const Eigen::Vector2f& uv, bool bilinear)
{
    return sample_color(tex, uv, bilinear).norm() / 255.f;
}

// 由法线 n 构造 TBN
// 思路：把世界坐标的“上”方向 ĵ=(0,1,0) 投影到切平面，残差归一化就是切线 t
// 当 n 平行于 ĵ 时投影长度 1-y² 为 0，此时随便取一组正交基兜底
static Eigen::Matrix3f build_tbn(const Eigen::Vector3f& n)
{
    const Eigen::Vector3f up(0.f, 1.f, 0.f);
    const float proj = up.dot(n);            // = n.y()
    const float len2 = 1.f - proj * proj;    // |ĵ - (ĵ·n)n|²

    Eigen::Vector3f t;
    if (len2 > 1e-8f)
        t = (up - proj * n) / std::sqrt(len2);
    else
        t = Eigen::Vector3f(1.f, 0.f, 0.f);  // n 与 ĵ 平行，切平面任意

    Eigen::Vector3f b = n.cross(t);
    b.normalize();

    Eigen::Matrix3f tbn;
    tbn.col(0) = t;
    tbn.col(1) = b;
    tbn.col(2) = n;
    return tbn;
}

// 高度扰动法线：在 uv 处做前向差分求高度梯度，得到切线空间法线 (-kh*du, -kh*dv, 1)
// 再乘 TBN 变回 view space。bump 和 displacement 共用这一段。
// strength = kh * kn：kh 控制“凹凸强度”，kn 是位移系数，这里作为小量一起乘在差分上
static Eigen::Vector3f perturb_normal(const Eigen::Vector3f& n, Texture* tex,
                                      const Eigen::Vector2f& uv, float strength)
{
    const Eigen::Matrix3f tbn = build_tbn(n);

    // 一个 texel 对应的 uv 步长
    const float du = tex ? 1.f / tex->width : 0.f;
    const float dv = tex ? 1.f / tex->height : 0.f;

    // 梯度用单通道 .x()，保持线性，避免三个通道做非线性组合
    const float h0 = sample_color(tex, uv, false).x();
    const float hu = sample_color(tex, Eigen::Vector2f(uv.x() + du, uv.y()), false).x();
    const float hv = sample_color(tex, Eigen::Vector2f(uv.x(), uv.y() + dv), false).x();

    const Eigen::Vector3f ln(-strength * (hu - h0), -strength * (hv - h0), 1.f);
    return (tbn * ln).normalized();
}

// 标准 Blinn-Phong：环境光 + 每盏点光源的漫反射与高光，返回 0~255 的颜色
static Eigen::Vector3f blinn_phong(const material& m, const Eigen::Vector3f& point,
                                   const Eigen::Vector3f& normal)
{
    // 环境光：与光源数量无关，只在循环外加一次
    Eigen::Vector3f color = m.ka.cwiseProduct(g_scene.amb);

    // 视线方向：从 shading point 指向相机
    Eigen::Vector3f view_dir = (g_scene.eye - point).normalized();

    for (const auto& l : g_scene.lights)
    {
        // 光照方向：从 shading point 指向光源
        Eigen::Vector3f light_dir = l.position - point;
        const float r2 = light_dir.squaredNorm();     // 距离衰减
        light_dir.normalize();

        // 漫反射
        Eigen::Vector3f diffuse = m.kd.cwiseProduct(l.intensity) * std::max(0.f, normal.dot(light_dir));

        // 高光：半程向量 h = normalize(L + V)
        Eigen::Vector3f half_dir = (view_dir + light_dir).normalized();
        Eigen::Vector3f specular = m.ks.cwiseProduct(l.intensity)
                                 * std::pow(std::max(0.f, normal.dot(half_dir)), m.p);

        color += (diffuse + specular) / r2;
    }

    return color * 255.f;
}

// ======================== 各个 fragment shader ========================
// 公共的照明都交给 blinn_phong，每个 shader 只写自己的差异部分

Eigen::Vector3f texture_fragment_shader(const fragment_shader_payload& payload)
{
    // 漫反射颜色来自纹理
    Eigen::Vector3f kd = sample_color(payload.texture, payload.tex_coords, payload.bilinear) / 255.f;
    material m{ Eigen::Vector3f(0.005f, 0.005f, 0.005f),  kd,
                Eigen::Vector3f(0.7937f, 0.7937f, 0.7937f), 150.f };
    return blinn_phong(m, payload.view_pos, payload.normal.normalized());
}

Eigen::Vector3f phong_fragment_shader(const fragment_shader_payload& payload)
{
    // 漫反射颜色来自顶点色（payload.color），不采样纹理
    material m{ Eigen::Vector3f(0.005f, 0.005f, 0.005f),  payload.color,
                Eigen::Vector3f(0.7937f, 0.7937f, 0.7937f), 150.f };
    return blinn_phong(m, payload.view_pos, payload.normal.normalized());
}

Eigen::Vector3f bump_fragment_shader(const fragment_shader_payload& payload)
{
    const float kh = 0.2f;
    const float kn = 0.1f;

    // 只扰动法线，不移动位置——这就是 bump 和 displacement 的唯一区别
    const Eigen::Vector3f normal = perturb_normal(payload.normal.normalized(), payload.texture,
        payload.tex_coords, kh * kn);

    // 漫反射颜色来自纹理（高度图同时当了 albedo，作业模板就是这么用的）
    Eigen::Vector3f kd = sample_color(payload.texture, payload.tex_coords, payload.bilinear) / 255.f;
    material m{ Eigen::Vector3f(0.005f, 0.005f, 0.005f),  kd,
                Eigen::Vector3f(0.7937f, 0.7937f, 0.7937f), 150.f };
    return blinn_phong(m, payload.view_pos, normal);
}

Eigen::Vector3f displacement_fragment_shader(const fragment_shader_payload& payload)
{
    // 与 bump 的区别只有一个：位移贴图把 shading point 沿法线推到新位置
    const float kh = 0.2f;
    const float kn = 0.1f;

    const Eigen::Vector3f normal = perturb_normal(payload.normal.normalized(), payload.texture,
                                                  payload.tex_coords, kh * kn);
    // 高度归一到 [0,1]，再乘 kn 作为位移量（单位与模型一致，模型高约 2.5）
    const float height = sample_height(payload.texture, payload.tex_coords, payload.bilinear);
    const Eigen::Vector3f point = payload.view_pos + kn * height * normal;

    Eigen::Vector3f kd = sample_color(payload.texture, payload.tex_coords, payload.bilinear) / 255.f;
    material m{ Eigen::Vector3f(0.005f, 0.005f, 0.005f),  kd,
                Eigen::Vector3f(0.7937f, 0.7937f, 0.7937f), 150.f };
    return blinn_phong(m, point, normal);
}

// ======================== 命令行 shader 表 ========================
// 加一个新 shader 只需要在这里加一行：关键字、函数、需要的纹理、是否双线性采样

struct shader_entry
{
    const char* name;
    std::function<Eigen::Vector3f(fragment_shader_payload)> fn;
    const char* texture;    // nullptr 表示沿用默认纹理
    bool bilinear;
};

static const std::vector<shader_entry> kShaders = {
    { "texture",                texture_fragment_shader,      "spot_texture.png",    false },
    { "texture_bilinear",       texture_fragment_shader,      "spot_texture.png",    true  },
    { "texture_small",          texture_fragment_shader,      "spot_texture_64.png", false },
    { "texture_small_bilinear", texture_fragment_shader,      "spot_texture_64.png", true  },
    { "normal",                 normal_fragment_shader,       nullptr,               false },
    { "phong",                  phong_fragment_shader,        nullptr,               false },
    { "bump",                   bump_fragment_shader,         "hmap.jpg",            false },
    { "displacement",           displacement_fragment_shader, "hmap.jpg",            false },
};

int main(int argc, const char** argv)
{
    std::vector<Triangle*> TriangleList;

    float angle = 180.0;
    bool command_line = false;

    std::string filename = "output.png";
    objl::Loader Loader;
    std::string obj_path = std::string(RASTERIZER_MODEL_DIR) + "/spot/";

    // Load .obj File
    bool loadout = Loader.LoadFile(obj_path + "spot_triangulated_good.obj");
    if (!loadout)
    {
        std::cerr << "Failed to load OBJ: " << obj_path << "spot_triangulated_good.obj\n";
        return 1;
    }
    for(auto mesh:Loader.LoadedMeshes)
    {
        for(int i=0;i<mesh.Vertices.size();i+=3)
        {
            Triangle* t = new Triangle();
            for(int j=0;j<3;j++)
            {
                t->setVertex(j,Vector4f(mesh.Vertices[i+j].Position.X,mesh.Vertices[i+j].Position.Y,mesh.Vertices[i+j].Position.Z,1.0));
                t->setNormal(j,Vector3f(mesh.Vertices[i+j].Normal.X,mesh.Vertices[i+j].Normal.Y,mesh.Vertices[i+j].Normal.Z));
                t->setTexCoord(j,Vector2f(mesh.Vertices[i+j].TextureCoordinate.X, mesh.Vertices[i+j].TextureCoordinate.Y));
            }
            TriangleList.push_back(t);
        }
    }

    rst::rasterizer r(700, 700);

    r.set_texture(Texture(obj_path + "spot_texture.png"));
    std::function<Eigen::Vector3f(fragment_shader_payload)> active_shader = displacement_fragment_shader;

    if (argc >= 2)
    {
        command_line = true;
        filename = std::string(argv[1]);

        if (argc == 3)
        {
            const std::string mode = argv[2];
            auto it = std::find_if(kShaders.begin(), kShaders.end(),
                                   [&](const shader_entry& e) { return mode == e.name; });
            if (it == kShaders.end())
            {
                std::cerr << "unknown shader: " << mode << "\n";
                return 1;
            }

            active_shader = it->fn;
            r.set_bilinear(it->bilinear);              // 采样方式随 payload 传给 shader
            if (it->texture != nullptr)                // normal/phong 沿用默认纹理
                r.set_texture(Texture(obj_path + it->texture));

            std::cout << "Rasterizing using the " << it->name << " shader\n";
        }
    }

    Eigen::Vector3f eye_pos = {0,0,10};

    r.set_vertex_shader(vertex_shader);
    r.set_fragment_shader(active_shader);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);
        cv::imwrite(filename, image);
        key = cv::waitKey(10);

        if (key == 'a' )
        {
            angle -= 10;
        }
        else if (key == 'd')
        {
            angle += 10;
        }

        angle -= 10;

    }
    return 0;
}
