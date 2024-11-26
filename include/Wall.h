#ifndef WALL_H
#define WALL_H

#include "Vector2D.h" // For position and normal vector handling
#include "Particle.h"

class Wall {
private:
    Vector2D start;   // Start point of the wall
    Vector2D end;     // End point of the wall
    Vector2D normal;  // Normal vector to the wall (used for reflections)

    // Helper function to calculate the normal vector
    void calculateNormal();

public:
    // Constructor
    Wall() = default;
    Wall(const Vector2D& start, const Vector2D& end);

    // Accessors
    [[nodiscard]] Vector2D getStart() const;
    [[nodiscard]] Vector2D getEnd() const;
    [[nodiscard]] Vector2D getNormal() const;

    // Collision detection: checks if a particle is colliding with the wall
    bool checkCollision(const Particle& particle, float dt) const;

    // Collision response: reflects particle velocity if a collision occurs
    void handleCollision(Particle& particle) const;

};

#endif // WALL_H