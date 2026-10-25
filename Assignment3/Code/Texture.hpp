//
// Created by LEI XU on 4/27/19.
//

#ifndef RASTERIZER_TEXTURE_H
#define RASTERIZER_TEXTURE_H
#include "global.hpp"
#include <algorithm>
#include <eigen3/Eigen/Eigen>
#include <opencv2/opencv.hpp>
#include <stdexcept>
class Texture{
private:
    cv::Mat image_data;

public:
    Texture(const std::string& name)
    {
        image_data = cv::imread(name);
        if (image_data.empty())
        {
            throw std::runtime_error("Failed to load texture: " + name);
        }
        cv::cvtColor(image_data, image_data, cv::COLOR_RGB2BGR);
        width = image_data.cols;
        height = image_data.rows;
    }

    int width, height;

    Eigen::Vector3f getColor(float u, float v)
    {
        // Clamp 寻址：越界的 UV 使用纹理边缘的颜色
        u = std::clamp(u, 0.0f, 1.0f);
        v = std::clamp(v, 0.0f, 1.0f);

        // UV 端点可能换算成 width/height，最大有效索引需要减一
        const int u_img = std::min(static_cast<int>(u * width), width - 1);
        const int v_img = std::min(static_cast<int>((1 - v) * height), height - 1);
        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        return Eigen::Vector3f(color[0], color[1], color[2]);
    }

    // 读取整数 texel 坐标（已保证在范围内）
    Eigen::Vector3f texel(int x, int y) const
    {
        auto color = image_data.at<cv::Vec3b>(y, x);
        return Eigen::Vector3f(color[0], color[1], color[2]);
    }

    // 双线性插值采样（Bonus）
    // 与 getColor 的唯一区别：不再直接截断取整，而是取"最近的 4 个 texel"按距离加权平均
    Eigen::Vector3f getColorBilinear(float u, float v)
    {
        u = std::clamp(u, 0.0f, 1.0f);
        v = std::clamp(v, 0.0f, 1.0f);

        // 1. UV -> 以 texel 为单位的浮点坐标（注意 v 翻转，且与 getColor 保持同一套约定）
        const float x = u * width;
        const float y = (1 - v) * height;

        // 2. 取左下角整点，并把小数部分当作插值权重
        const int x0 = static_cast<int>(x);
        const int y0 = static_cast<int>(y);
        const float s = x - x0;   // 水平方向权重，∈[0,1)
        const float t = y - y0;   // 垂直方向权重，∈[0,1)

        // 3. 四个邻居，边界要 clamp，否则 u=1/v=0 会越界
        const int x1 = std::min(x0 + 1, width - 1);
        const int y1 = std::min(y0 + 1, height - 1);
        const int xc = std::min(x0, width - 1);
        const int yc = std::min(y0, height - 1);

        const Eigen::Vector3f c00 = texel(xc, yc);
        const Eigen::Vector3f c10 = texel(x1, yc);
        const Eigen::Vector3f c01 = texel(xc, y1);
        const Eigen::Vector3f c11 = texel(x1, y1);

        // 4. 先水平插值两行，再垂直插值一次
        const Eigen::Vector3f top = c00 * (1 - s) + c10 * s;
        const Eigen::Vector3f bottom = c01 * (1 - s) + c11 * s;
        return top * (1 - t) + bottom * t;
    }

};
#endif //RASTERIZER_TEXTURE_H
