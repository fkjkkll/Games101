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

};
#endif //RASTERIZER_TEXTURE_H
