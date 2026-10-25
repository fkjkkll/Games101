#include <iostream>
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
    Eigen::Vector3f return_color = (payload.normal.head<3>() + Eigen::Vector3f(1.0f, 1.0f, 1.0f)) / 2.f;
    Eigen::Vector3f result;
    result << return_color.x() * 255, return_color.y() * 255, return_color.z() * 255;
    return result;
}

static Eigen::Vector3f reflect(const Eigen::Vector3f& vec, const Eigen::Vector3f& axis)
{
    auto costheta = vec.dot(axis);
    return (2 * costheta * axis - vec).normalized();
}

struct light
{
    Eigen::Vector3f position;
    Eigen::Vector3f intensity;
};

Eigen::Vector3f texture_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f return_color = {0, 0, 0};
    if (payload.texture)
    {
        auto tex_coord = payload.tex_coords;
        return_color = payload.texture->getColor(tex_coord.x(), tex_coord.y());
        // TODO: Get the texture value at the texture coordinates of the current fragment
    }
    Eigen::Vector3f texture_color;
    texture_color << return_color.x(), return_color.y(), return_color.z();

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = texture_color / 255.f;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = texture_color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};

    Vector3f viewVector = (eye_pos - point).normalized();
    for (auto& light : lights)
    {
        auto lightVector = light.position - point;
        auto lightUnitVector = lightVector.normalized();
        auto distSqr = lightVector.squaredNorm();

        // 漫反射
        auto reflect = MAX(normal.dot(lightUnitVector), 0);
        Vector3f diffuse = kd.cwiseProduct(light.intensity / distSqr) * reflect;

        // 高光项
        auto halfUnitVector = (lightUnitVector + viewVector).normalized();
        auto mirror = MAX(normal.dot(halfUnitVector), 0);
        auto specular = ks.cwiseProduct(light.intensity / distSqr) * pow(mirror, p);

        // 环境光
        auto ambient = ka.cwiseProduct(amb_light_intensity);

        result_color += diffuse;
        result_color += specular;
        result_color += ambient;

        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
    }

    return result_color * 255.f;
}

Eigen::Vector3f phong_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};

    Vector3f viewVector = (eye_pos - point).normalized();
    for (auto& light : lights)
    {
        auto lightVector = light.position - point;
        auto lightUnitVector = lightVector.normalized();
        auto distSqr = lightVector.squaredNorm();

        // 漫反射
        auto reflect = MAX(normal.dot(lightUnitVector), 0);
        Vector3f diffuse = kd.cwiseProduct(light.intensity / distSqr) * reflect;

        // 高光项
        auto halfUnitVector = (lightUnitVector + viewVector).normalized();
        auto mirror = MAX(normal.dot(halfUnitVector), 0);
        auto specular = ks.cwiseProduct(light.intensity / distSqr) * pow(mirror, p);

        // 环境光
        auto ambient = ka.cwiseProduct(amb_light_intensity);

        result_color += diffuse;
        result_color += specular;
        result_color += ambient;

        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
    }

    return result_color * 255.f;
}

Eigen::Vector3f displacement_fragment_shader(const fragment_shader_payload& payload)
{
    
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;

    // 0. 纹理颜色当作 kd：位移的细节主要靠 albedo 的变化体现，不能用 payload.color
    Eigen::Vector3f texture_color(0, 0, 0);
    if (payload.texture)
        texture_color = payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()) / 255.f;
    Eigen::Vector3f kd = texture_color;

    // 1. 构造切线空间 TBN（与 bump_fragment_shader 完全相同）
    normal = normal.normalized();
    const Vector3f up(0, 1, 0);
    const float proj = normal.dot(up); // = n.y()
    const float s2 = 1.0f - proj * proj; // 1 - y² = x² + z²

    Vector3f t;
    if (s2 > 1e-6f)
        t = (up - proj * normal) / sqrt(s2); // 别忘了归一化
    else
        t = Vector3f(1, 0, 0); // n 平行 ±Y 时的退化兜底

    const Vector3f b = normal.cross(t); // 副切线
    Matrix3f tbn;
    tbn << t, b, normal;

    // 2. 高度图：差分 → 扰动法线；高度值 → 位置位移
    if (payload.texture)
    {
        Texture* const tex = payload.texture;
        const float deltaU = 1.0f / static_cast<float>(tex->width);
        const float deltaV = 1.0f / static_cast<float>(tex->height);
        const Vector2f uv = payload.tex_coords;

        const float h0 = tex->getColor(uv.x(), uv.y()).x();
        const float hu = tex->getColor(uv.x() + deltaU, uv.y()).x();
        const float hv = tex->getColor(uv.x(), uv.y() + deltaV).x();

        const float du = kh * kn * (hu - h0); // 差分是小量，所以乘 kh*kn
        const float dv = kh * kn * (hv - h0);

        const Vector3f ln(-du, -dv, 1.0f);
        normal = (tbn * ln).normalized(); // 这一行和 bump 一模一样

        // 位移量 = kn * h。量纲是这题最大的坑：getColor 返回 0~255，
        // 直接写 kn * h 会挤出 0~25 个单位，而模型本身只有 2.5 大，整只牛会被炸开。
        // 所以先把 h 归一化到 0~1（想要更夸张的鼓起，把 kn 调到 0.5~1.0 即可）。
        const float height = tex->getColor(uv.x(), uv.y()).norm() / 255.f;
        point += kn * height * normal; // 注意：沿"扰动后"的法线挤出
    }

    // 3. Blinn-Phong 光照（公式与 phong_fragment_shader 一致，但 point 已经移动过）
    Eigen::Vector3f result_color = {0, 0, 0};
    const Vector3f viewVector = (eye_pos - point).normalized();

    for (auto& light : lights)
    {
        auto lightVector = light.position - point;
        auto lightUnitVector = lightVector.normalized();
        auto distSqr = lightVector.squaredNorm();

        // 漫反射
        auto reflect = MAX(normal.dot(lightUnitVector), 0);
        Vector3f diffuse = kd.cwiseProduct(light.intensity / distSqr) * reflect;

        // 高光项
        auto halfUnitVector = (lightUnitVector + viewVector).normalized();
        auto mirror = MAX(normal.dot(halfUnitVector), 0);
        auto specular = ks.cwiseProduct(light.intensity / distSqr) * pow(mirror, p);

        // 环境光
        auto ambient = ka.cwiseProduct(amb_light_intensity);

        result_color += diffuse;
        result_color += specular;
        result_color += ambient;
    }

    return result_color * 255.f;
}

Eigen::Vector3f bump_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f return_color = { 0, 0, 0 };
    if (payload.texture)
    {
        auto tex_coord = payload.tex_coords;
        return_color = payload.texture->getColor(tex_coord.x(), tex_coord.y());
    }
    Eigen::Vector3f texture_color;
    texture_color << return_color.x(), return_color.y(), return_color.z();

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = texture_color / 255.f;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    // bump
    float kh = 0.2, kn = 0.1;

    // 1. 构造切线空间 TBN
    // 需要 t、b、n 三个互相正交的"单位"向量：t 或 b 的长度不为 1 时，
    // TBN 会把 ln 里的扰动一起缩放，凹凸强度就会随法线朝向变化。
    normal = normal.normalized(); // 下面假设 |n| = 1，保险起见再归一化一次

    const Vector3f up(0, 1, 0); // 参考的"向上"轴
    const float proj = normal.dot(up); // = n.y()
    const float s2 = 1.0f - proj * proj; // = |up - proj*n|² = x² + z²

    Vector3f t;
    if (s2 > 1e-6f)
    {
        t = (up - proj * normal) / sqrt(s2); // 归一化：除以 |up - proj*n| = √s² = s
    }
    else
    {
        t = Vector3f(1, 0, 0); // n 平行 ±Y 轴时退化，任选一个垂直方向兜底
    }

    const Vector3f b = normal.cross(t); // n、t 均为单位正交 ⇒ b 天然是单位向量
    Matrix3f tbn;
    tbn << t, b, normal;

    // 2. 用纹理的 R 通道当作高度场，前向差分求梯度
    if (payload.texture)
    {
        Texture* const tex = payload.texture;
        const float deltaU = 1.0f / static_cast<float>(tex->width); // 一个 texel 的 UV 步长
        const float deltaV = 1.0f / static_cast<float>(tex->height);
        const Vector2f uv = payload.tex_coords;

        const float h0 = tex->getColor(uv.x(), uv.y()).x();
        const float hu = tex->getColor(uv.x() + deltaU, uv.y()).x();
        const float hv = tex->getColor(uv.x(), uv.y() + deltaV).x();

        const float du = kh * kn * (hu - h0);
        const float dv = kh * kn * (hv - h0);

        const Vector3f ln(-du, -dv, 1.0f);
        normal = (tbn * ln).normalized();
    }

    // TODO: Implement bump mapping here
    // Let n = normal = (x, y, z)
    // Vector t = (x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z))
    // Vector b = n cross product t
    // Matrix TBN = [t b n]
    // dU = kh * kn * (h(u+1/w,v)-h(u,v))
    // dV = kh * kn * (h(u,v+1/h)-h(u,v))
    // Vector ln = (-dU, -dV, 1)
    // Normal n = normalize(TBN * ln)

    Eigen::Vector3f result_color = { 0, 0, 0 };

    Vector3f viewVector = (eye_pos - point).normalized();
    for (auto& light : lights)
    {
        auto lightVector = light.position - point;
        auto lightUnitVector = lightVector.normalized();
        auto distSqr = lightVector.squaredNorm();

        // 漫反射
        auto reflect = MAX(normal.dot(lightUnitVector), 0);
        Vector3f diffuse = kd.cwiseProduct(light.intensity / distSqr) * reflect;

        // 高光项
        auto halfUnitVector = (lightUnitVector + viewVector).normalized();
        auto mirror = MAX(normal.dot(halfUnitVector), 0);
        auto specular = ks.cwiseProduct(light.intensity / distSqr) * pow(mirror, p);

        // 环境光
        auto ambient = ka.cwiseProduct(amb_light_intensity);

        result_color += diffuse;
        result_color += specular;
        result_color += ambient;
    }

    return result_color * 255.f;
}

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

    auto texture_path = "spot_texture.png";
    //auto texture_path = "hmap.jpg";
    r.set_texture(Texture(obj_path + texture_path));

    std::function<Eigen::Vector3f(fragment_shader_payload)> active_shader = displacement_fragment_shader;

    if (argc >= 2)
    {
        command_line = true;
        filename = std::string(argv[1]);

        if (argc == 3 && std::string(argv[2]) == "texture")
        {
            std::cout << "Rasterizing using the texture shader\n";
            active_shader = texture_fragment_shader;
            texture_path = "spot_texture.png";
            r.set_texture(Texture(obj_path + texture_path));
        }
        else if (argc == 3 && std::string(argv[2]) == "normal")
        {
            std::cout << "Rasterizing using the normal shader\n";
            active_shader = normal_fragment_shader;
        }
        else if (argc == 3 && std::string(argv[2]) == "phong")
        {
            std::cout << "Rasterizing using the phong shader\n";
            active_shader = phong_fragment_shader;
        }
        else if (argc == 3 && std::string(argv[2]) == "bump")
        {
            std::cout << "Rasterizing using the bump shader\n";
            active_shader = bump_fragment_shader;
            texture_path = "hmap.jpg";   // bump/displacement 要用高度图
            r.set_texture(Texture(obj_path + texture_path));
        }
        else if (argc == 3 && std::string(argv[2]) == "displacement")
        {
            std::cout << "Rasterizing using the displacement shader\n";
            active_shader = displacement_fragment_shader;
            texture_path = "hmap.jpg";   // bump/displacement 要用高度图
            r.set_texture(Texture(obj_path + texture_path));
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
