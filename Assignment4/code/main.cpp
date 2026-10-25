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

// ==================== 抗锯齿 ====================

enum class AntiAliasMode
{
    None,        // 无抗锯齿：把浮点坐标直接截断成整数像素
    Splat,       // 双线性 splat：把一个采样点的能量分摊给周围 4 个像素
    Supersample, // 3x3 超采样：每个像素取 9 个子采样点，统计被覆盖的比例
};

// 想换一种画法，只改这一行
constexpr AntiAliasMode kAntiAliasMode = AntiAliasMode::Splat;

/// <summary>
/// 点到线段的最短距离（超采样判断"子采样点是否被曲线覆盖"时用）
/// （注意，这里不是点到直线的距离，如果去掉clamp就是了）
/// </summary>
static float distance_to_segment(const cv::Point2f &p, const cv::Point2f &a, const cv::Point2f &b)
{
    const cv::Point2f ab = b - a;
    const cv::Point2f ap = p - a;
    const float denom = ab.dot(ab);
    float t = (denom > 1e-12f) ? ap.dot(ab) / denom : 0.0f;
    t = std::min<float>(1.0f, std::max<float>(0.0f, t)); // 投影点夹在线段内
    return cv::norm(p - (a + ab * t));
}

// 把曲线离散成折线：三种画法共用同一批采样点，对照才公平
static std::vector<cv::Point2f> sample_curve(const std::vector<cv::Point2f> &control_points)
{
    std::vector<cv::Point2f> polyline;
    for (double t = 0.0; t <= 1.0; t += 0.001)
        polyline.push_back(recursive_bezier(control_points, t));
    return polyline;
}

// 把覆盖率缓冲合成到图像：alpha = min(1, coverage)，与黑底混合（写绿通道）
static void composite_coverage(const std::vector<float> &coverage, cv::Mat &window)
{
    for (int y = 0; y < window.rows; ++y) {
        for (int x = 0; x < window.cols; ++x) {
            const float alpha = std::min<float>(1.0f, coverage[static_cast<size_t>(y) * window.cols + x]);
            if (alpha <= 0.0f) continue;

            auto& pixel = window.at<cv::Vec3b>(y, x);
            // 用显式模板参数，避免 windows.h 的 max 宏干扰
            pixel[1] = std::max<uchar>(pixel[1], static_cast<uchar>(255.0f * alpha + 0.5f)); // 省事版四舍五入
        }
    }
}

// 画法 1：无抗锯齿 —— 一个采样点只影响一个像素，非 0 即 255
static void rasterize_no_aa(const std::vector<cv::Point2f> &polyline, cv::Mat &window)
{
    for (const auto &point : polyline) {
        const int x = static_cast<int>(point.x);
        const int y = static_cast<int>(point.y);
        if (x < 0 || x >= window.cols || y < 0 || y >= window.rows) continue;
        window.at<cv::Vec3b>(y, x)[1] = 255;
    }
}

// 画法 2：双线性 splat —— 覆盖率是连续量，本质是积分，积分就要累加
static void rasterize_splat(const std::vector<cv::Point2f> &polyline, cv::Mat &window)
{
    const int width = window.cols;
    const int height = window.rows;

    // ☆ 抗锯齿的本质是积分，积分就要累加 ☆
    // 覆盖率缓冲：每个像素被曲线"覆盖"了多少，取值 0~1 之间连续变化
    std::vector<float> coverage(static_cast<size_t>(width) * height, 0.0f);

    for (const auto &point : polyline) {
        // 拆成整数格点 + 小数部分，小数部分就是点在像素内的相对位置
        const int x0 = static_cast<int>(std::floor(point.x));
        const int y0 = static_cast<int>(std::floor(point.y));
        const float fx = point.x - static_cast<float>(x0);
        const float fy = point.y - static_cast<float>(y0);

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

    composite_coverage(coverage, window);
}

// 画法 3：3x3 覆盖率采样（MSAA 语义；纯色时结果与 SSAA 等价）—— 每个像素取 9 个点，判断子采样点是否落在曲线（半宽 0.5）内
static void rasterize_supersample(const std::vector<cv::Point2f> &polyline, cv::Mat &window)
{
    const int width = window.cols;
    const int height = window.rows;
    constexpr int n = 3;                                           // 3x3 子采样
    constexpr float offsets[n] = {1.0f / 6.0f, 0.5f, 5.0f / 6.0f}; // 子采样点在像素内的位置
    constexpr float half_width = 0.5f;                             // 线宽 1 像素 => 半宽 0.5

    // 每个像素用 9 个 bit 记录哪些子采样点被覆盖，多条线段取并集
    std::vector<unsigned short> hits(static_cast<size_t>(width) * height, 0);

    for (size_t i = 1; i < polyline.size(); ++i) {
        const cv::Point2f a = polyline[i - 1];
        const cv::Point2f b = polyline[i];

        // 只处理这条线段附近的像素
        // >> 线段与半径0.5的圆形做闵可夫斯基和，因此是实际线段AABB再向外扩0.5，向上取整就要多扩一圈像素 <<
        // 严格的推导是（以x举例）：x + 5/6 >= minx - 0.5 => x >= minx - 4/3 = ceil(minx - 4/3) 我们的宽松界是 floor(minx) - 1
        //                      x + 1/6 <= maxx + 0.5 => x <= maxx + 1/3 = floor(maxx + 1/3) 我们的宽松界是 ceil(maxx) + 1
        const int xmin = std::max(0, static_cast<int>(std::floor(std::min(a.x, b.x))) - 1);
        const int xmax = std::min(width - 1, static_cast<int>(std::ceil(std::max(a.x, b.x))) + 1);
        const int ymin = std::max(0, static_cast<int>(std::floor(std::min(a.y, b.y))) - 1);
        const int ymax = std::min(height - 1, static_cast<int>(std::ceil(std::max(a.y, b.y))) + 1);

        for (int y = ymin; y <= ymax; ++y) {
            for (int x = xmin; x <= xmax; ++x) {
                unsigned short mask = 0;
                for (int iy = 0; iy < n; ++iy)
                    for (int ix = 0; ix < n; ++ix)
                    {
                        const cv::Point2f sample(x + offsets[ix], y + offsets[iy]);
                        if (distance_to_segment(sample, a, b) <= half_width) // 覆盖测试
                            mask |= static_cast<unsigned short>(1u << (iy * n + ix));
                    }
                hits[static_cast<size_t>(y) * width + x] |= mask;
            }
        }
    }

    std::vector<float> coverage(static_cast<size_t>(width) * height, 0.0f);
    for (size_t i = 0; i < coverage.size(); ++i) {
        int covered = 0;
        for (unsigned short m = hits[i]; m; m >>= 1) covered += (m & 1u);
        coverage[i] = static_cast<float>(covered) / (n * n);
    }

    // coverage 里的值其实不需要 clamp，里面的取值是固定的：0, 1/9, 2/9, ..., 1
    composite_coverage(coverage, window);
}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window, AntiAliasMode mode)
{
    if (control_points.empty()) return;

    const std::vector<cv::Point2f> polyline = sample_curve(control_points);

    switch (mode) {
        case AntiAliasMode::None:
            rasterize_no_aa(polyline, window);
            break;
        case AntiAliasMode::Supersample:
            rasterize_supersample(polyline, window);
            break;
        case AntiAliasMode::Splat:
        default:
            rasterize_splat(polyline, window);
            break;
    }
}

// 作业要求的入口：用常量 kAntiAliasMode 选择画法
void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window)
{
    bezier(control_points, window, kAntiAliasMode);
}

// 左上角显示当前抗锯齿模式，底部显示快捷键
static void draw_hud(cv::Mat &image, AntiAliasMode mode)
{
    const char *name =
        (mode == AntiAliasMode::None) 
        ? "AA: None" 
        : (mode == AntiAliasMode::Splat)
            ? "AA: Splat" 
            : "AA: Supersample 3x3";
    cv::putText(image, name, {12, 30}, cv::FONT_HERSHEY_SIMPLEX, 0.7, {200, 200, 200}, 1, cv::LINE_AA);
    cv::putText(image, "1=None  2=Splat  3=Supersample   r=reset   s=save   Esc=quit", {12, 688},
                cv::FONT_HERSHEY_SIMPLEX, 0.45, {110, 110, 110}, 1, cv::LINE_AA);
}

int main() 
{
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);
    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    AntiAliasMode mode = kAntiAliasMode;
    int key = -1;
    int last_count = 0;
    bool need_save = false;

    while (key != 27) {
        key = cv::waitKey(20);

        // 快捷键：1/2/3 切换抗锯齿模式，r 清空控制点，s 保存当前画面
        if (key == '1') { mode = AntiAliasMode::None; need_save = true; }
        else if (key == '2') { mode = AntiAliasMode::Splat; need_save = true; }
        else if (key == '3') { mode = AntiAliasMode::Supersample; need_save = true; }
        else if (key == 'r' || key == 'R') { control_points.clear(); last_count = 0; }
        else if (key == 's' || key == 'S') { need_save = true; }

        // 每帧从空白画布重画：切换模式后立刻能看到区别，也不会把旧图叠加进来
        cv::Mat canvas(700, 700, CV_8UC3, cv::Scalar(0));
        // BGR2RGB: 红绿蓝 -> 蓝绿红
        cv::cvtColor(canvas, canvas, cv::COLOR_BGR2RGB);

        // 控制多边形：把相邻控制点连起来，可以直观看到曲线被"拉向"控制点的效果
        for (size_t i = 1; i < control_points.size(); ++i)
            cv::line(canvas, control_points[i - 1], control_points[i], {140, 140, 140}, 1, cv::LINE_AA);

        for (auto &point : control_points)
            cv::circle(canvas, point, 3, {255, 255, 255}, 3);

        if (control_points.size() == 4) {
            //naive_bezier(control_points, canvas);   // 想看"无抗锯齿"对照可以打开这行（画红通道）
            bezier(control_points, canvas, mode);
        }

        // 点满 4 个点、或刚切换过模式时，把当前画面存盘
        if (control_points.size() == 4 && (last_count < 4 || need_save)) {
            cv::imwrite("my_bezier_curve.png", canvas);
            need_save = false;
        }
        last_count = static_cast<int>(control_points.size());

        draw_hud(canvas, mode);
        cv::imshow("Bezier Curve", canvas);
    }

    return 0;
}
