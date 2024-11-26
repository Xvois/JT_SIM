//
// Created by Sonny Parker on 25/11/2024.
//

#ifndef QAUD_H
#define QAUD_H


class Quad {
public:
    double x, y, width, height;

    Quad(double x, double y, double width, double height)
        : x(x), y(y), width(width), height(height) {}

    [[nodiscard]] bool contains(double px, double py) const {
        return (px >= x - width / 2 && px <= x + width / 2 &&
                py >= y - height / 2 && py <= y + height / 2);
    }

    [[nodiscard]] bool intersects(const Quad& range) const {
        return !(range.x - range.width / 2 > x + width / 2 ||
                 range.x + range.width / 2 < x - width / 2 ||
                 range.y - range.height / 2 > y + height / 2 ||
                 range.y + range.height / 2 < y - height / 2);
    }

    [[nodiscard]] double getVolume() const {
        return width * height;
    }
};

#endif //QAUD_H
