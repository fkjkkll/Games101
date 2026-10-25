#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 4) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, double t)
{
    size_t count = control_points.size();
    if (count == 1)
        return control_points[0];

    std::vector<cv::Point2f> next_control_points;
    next_control_points.reserve(count - 1);

    for (size_t i = 1; i < count; ++i) {
        auto new_point = control_points[i - 1] * (1 - t) + control_points[i] * t;
        next_control_points.push_back(new_point);
    }

    return recursive_bezier(next_control_points, t);
}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    if (control_points.empty()) return;

    const int width = window.cols;
    const int height = window.rows;

    // ☆ 抗锯齿的本质是积分，积分就要累加 ☆
    // 覆盖率缓冲：每个像素被曲线"覆盖"了多少，取值 0~1 之间连续变化
    std::vector<float> coverage(static_cast<size_t>(width) * height, 0.0f);

    for (double t = 0.0; t <= 1.0; t += 0.001) {
        cv::Point2f point = recursive_bezier(control_points, t);

        // 拆成整数格点 + 小数部分，小数部分就是点在像素内的相对位置
        int x0 = static_cast<int>(std::floor(point.x));
        int y0 = static_cast<int>(std::floor(point.y));
        float fx = point.x - static_cast<float>(x0);
        float fy = point.y - static_cast<float>(y0);

        // 双线性权重：把这一点的能量按距离分摊给周围 4 个像素
        const float w00 = (1.0f - fx) * (1.0f - fy);
        const float w10 = fx * (1.0f - fy);
        const float w01 = (1.0f - fx) * fy;
        const float w11 = fx * fy;

        if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height)
            coverage[static_cast<size_t>(y0) * width + x0] += w00;
        if (x0 + 1 >= 0 && x0 + 1 < width && y0 >= 0 && y0 < height)
            coverage[static_cast<size_t>(y0) * width + x0 + 1] += w10;
        if (x0 >= 0 && x0 < width && y0 + 1 >= 0 && y0 + 1 < height)
            coverage[static_cast<size_t>(y0 + 1) * width + x0] += w01;
        if (x0 + 1 >= 0 && x0 + 1 < width && y0 + 1 >= 0 && y0 + 1 < height)
            coverage[static_cast<size_t>(y0 + 1) * width + x0 + 1] += w11;
    }

    // 合成：alpha = min(1, coverage)，与黑底做线性混合
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float alpha = std::min<float>(1.0f, coverage[static_cast<size_t>(y) * width + x]);
            if (alpha <= 0.0f) continue;

            auto &pixel = window.at<cv::Vec3b>(y, x);
            uchar value = static_cast<uchar>(255.0f * alpha + 0.5f);
            // 用显式模板参数，避免 windows.h 的 max 宏干扰
            pixel[1] = std::max<uchar>(pixel[1], value);
        }
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    // BGR2RGB: 红绿蓝 -> 蓝绿红
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        // 控制多边形：把相邻控制点连起来，可以直观看到曲线被"拉向"控制点的效果
        for (size_t i = 1; i < control_points.size(); ++i)
        {
            cv::line(window, control_points[i - 1], control_points[i], {140, 140, 140}, 1, cv::LINE_AA);
        }

        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == 4) 
        {
            //naive_bezier(control_points, window);
            bezier(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            key = cv::waitKey(0);

            return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}
