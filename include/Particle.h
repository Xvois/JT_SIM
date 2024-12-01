//
// Created by Sonny Parker on 21/11/2024.
//
#ifndef PARTICLE_H
#define PARTICLE_H
#include "Vector2D.h"

class Particle
{
    static int next_id; // Static member to keep track of the next available id

protected:
    int id;
    Vector2D position;
    Vector2D velocity;
    Vector2D acceleration;
    float mass;

public:
    Particle() : id(next_id++), position(), velocity(), mass(0)
    {
    }

    Particle(const Vector2D& position, const Vector2D& velocity, float mass)
        : id(next_id++), position(position), velocity(velocity), mass(mass)
    {
    }


    virtual ~Particle() = default;
    virtual void update(float dt);
    void applyImpulse(const Vector2D& impulse);
    [[nodiscard]] Vector2D predictPosition(float dt) const;
    [[nodiscard]] Vector2D getPosition() const
    {
        return position;
    }
    void setPosition(const Vector2D& position)
    {
        this->position = position;
    }
    [[nodiscard]] Vector2D getVelocity() const
    {
        return velocity;
    }
    [[nodiscard]] Vector2D getAcceleration() const {
        return acceleration;
    }
    [[nodiscard]] float getMass() const
    {
        return mass;
    }
    [[nodiscard]] float getKineticEnergy() const
    {
        return 0.5f * mass * Vector2D::magnitudeSquared(velocity);
    }
    [[nodiscard]] int getId() const
    {
        return id;
    }

    // Needed for use in std::vector
    bool operator==(const Particle& other) const
    {
        return this->id == other.id;
    }
};

#endif //PARTICLE_H
