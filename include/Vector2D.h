//
// Created by Sonny Parker on 21/11/2024.
//

#ifndef VECTOR2D_H
#define VECTOR2D_H
#include <cmath>

class Vector2D
{
    public:
    double x;
    double y;

    // Default constructor
    Vector2D() : x(0), y(0) {}

    // Constructor with x and y components
    Vector2D(const double x, const double y) : x(x), y(y) {}

    // Copy constructor
    Vector2D(const Vector2D& other) : x(other.x), y(other.y) {}


    Vector2D operator+(Vector2D const& target) const
    {
        return { this->x + target.x, this->y + target.y };
    }

    Vector2D operator-(Vector2D const& target) const
    {
        return { this->x - target.x, this->y - target.y };
    }

    Vector2D operator-() const {
        return {-x, -y};
    }

    Vector2D operator*(const double scalar) const {
        return {x * scalar, y * scalar};
    }

    friend Vector2D operator*(const double scalar, const Vector2D& vec) {
        return {vec.x * scalar, vec.y * scalar};
    }

    Vector2D operator/(const double scalar) const
    {
        return {x / scalar, y / scalar};
    }

    Vector2D operator +=(const Vector2D& target)
    {
        this->x += target.x;
        this->y += target.y;
        return *this;
    }

    Vector2D operator -=(const Vector2D& target)
    {
        this->x -= target.x;
        this->y -= target.y;
        return *this;
    }

    static float magnitudeSquared(const Vector2D& velocity)
    {
        return velocity.x * velocity.x + velocity.y * velocity.y;
    }


    static double dot(Vector2D const& v1, Vector2D const& v2)
    {
        return v1.x * v2.x + v1.y * v2.y;
    }

    static Vector2D cross(Vector2D const& v1, Vector2D const& v2)
    {
        return {v1.x * v2.y, -v1.y * v2.x};
    }

    static double magnitude(Vector2D const& v1)
    {
        return sqrt(v1.x * v1.x + v1.y * v1.y);
    }

    static Vector2D normalize(Vector2D const& v1)
    {
        const double mag = magnitude(v1);
        return {v1.x / mag, v1.y / mag};
    }

};

#endif //VECTOR2D_H
