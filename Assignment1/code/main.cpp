#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

	auto rotationRadian = rotation_angle / 180.0 * MY_PI;
    Eigen::Matrix4f rotationMatrix;
    rotationMatrix << 
        cos(rotationRadian), -sin(rotationRadian), 0, 0,
        sin(rotationRadian), cos(rotationRadian), 0, 0,
        0, 0, 1, 0,
		0, 0, 0, 1;
	model = rotationMatrix * model;

    return model;
}

/// <summary>
/// 罗德里格斯旋转公式
/// 参考1: https://www.bilibili.com/video/BV1X7411F744?spm_id_from=333.788.videopod.episodes&vd_source=5def309c50dba10ca2a52c71628df16d&p=4
/// 参考2: https://zhuanlan.zhihu.com/p/451579313
/// 参考3: https://immortalqx.github.io/2021/07/08/slam-learning-2/
/// 参考4: https://www.bilibili.com/video/BV1qV411X79j/?spm_id_from=333.337.search-card.all.click&vd_source=5def309c50dba10ca2a52c71628df16d
/// </summary>
Eigen::Matrix4f get_rotation(Vector3f axis, float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    auto rotationRadian = rotation_angle / 180.0 * MY_PI;
    Eigen::Matrix4f rotationMatrix;
    rotationMatrix <<
        cos(rotationRadian) + axis.x() * axis.x() * (1 - cos(rotationRadian)),
        axis.x()* axis.y()* (1 - cos(rotationRadian)) - axis.z() * sin(rotationRadian),
        axis.x()* axis.z()* (1 - cos(rotationRadian)) + axis.y() * sin(rotationRadian),
        0,
        axis.y()* axis.x()* (1 - cos(rotationRadian)) + axis.z() * sin(rotationRadian),
        cos(rotationRadian) + axis.y() * axis.y() * (1 - cos(rotationRadian)),
        axis.y()* axis.z()* (1 - cos(rotationRadian)) - axis.x() * sin(rotationRadian),
        0,
        axis.z()* axis.x()* (1 - cos(rotationRadian)) - axis.y() * sin(rotationRadian),
        axis.z()* axis.y()* (1 - cos(rotationRadian)) + axis.x() * sin(rotationRadian),
        cos(rotationRadian) + axis.z() * axis.z() * (1 - cos(rotationRadian)),
        0,
        0, 0, 0, 1;
    model = rotationMatrix * model;
    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    double t = tan(eye_fov / 2.0 / 180.0 * MY_PI) * abs(zNear);
    double b = -t;
	double r = aspect_ratio * t;
    double l = -r;

	// 近平面和远平面都应该是负数
    zNear = -abs(zNear);
    zFar = -abs(zFar);

    Eigen::Matrix4f transitionMatrix;
	transitionMatrix << 
        1, 0, 0, -(r + l) / 2.0,
        0, 1, 0, -(t + b) / 2.0,
        0, 0, 1, -(zNear + zFar) / 2.0,
        0, 0, 0, 1;
	Eigen::Matrix4f scaleMatrix;
    scaleMatrix << 
        2.0 / (r - l), 0, 0, 0,
        0, 2.0 / (t - b), 0, 0,
        0, 0, 2.0 / (zNear - zFar), 0,
		0, 0, 0, 1;
    Eigen::Matrix4f persToOrthMatrix;
    persToOrthMatrix << 
        zNear, 0, 0, 0,
        0, zNear, 0, 0,
        0, 0, zNear + zFar, -zNear * zFar,
		0, 0, 1, 0;
	projection = scaleMatrix * transitionMatrix * persToOrthMatrix * projection;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    float axisAngle = 0;
    Eigen::Vector3f axis = { -1, 0, 0 };
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 10};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_rotation(axis, axisAngle) * get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_rotation(axis, axisAngle) * get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }

        if (key == 'w') {
            axisAngle += 10;
        }
        else if (key == 's') {
            axisAngle -= 10;
        }
    }

    return 0;
}
