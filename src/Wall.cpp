//
// Created by Sonny Parker on 21/11/2024.
//

#include "../include/Wall.h"

#include <iostream>

void Wall::calculateNormal()
{
    // Stores the normal of the current wall
    Vector2D delta = this-> end - this-> start;
    this-> normal = Vector2D::normalize(Vector2D(-delta.y, delta.x));
}


Wall::Wall(const Vector2D& start, const Vector2D& end)
    : start(start), end(end)
{
    this-> start = start;
    this-> end = end;
    this->calculateNormal();
}

Vector2D Wall::getStart() const
{
    return this-> start;
}

Vector2D Wall::getEnd() const
{
    return this-> end;
}

Vector2D Wall::getNormal() const
{
    return this-> normal;
}

/**
 * @param particle Particle to check for collision
 * @param dt Timestep over which to check
 * @return Boolean indicating if a collision will occur
 */
bool Wall::checkCollision(const Particle& particle, float dt) const
{
    Vector2D pos = particle.getPosition();               // Current position
    Vector2D predicted_pos = particle.predictPosition(dt); // Predicted position

    // Calculate the vector from the current position to the predicted position
    Vector2D movement = predicted_pos - pos;

    // Calculate the normal of the wall
    Vector2D wall_vec = this->end - this->start;
    Vector2D wall_normal = Vector2D::normalize(Vector2D(-wall_vec.y, wall_vec.x));

    // Calculate the distance from the particle to the wall at the current position
    float curr_dist = Vector2D::dot(pos - this->start, wall_normal);

    // Calculate the distance from the particle to the wall at the predicted position
    float predict_dist = Vector2D::dot(predicted_pos - this->start, wall_normal);

    // Check if there is a sign change in the distance (collision with infinite wall)
    if (curr_dist * predict_dist >= 0)
        return false;

    // Calculate the intersection point of the particle's path with the wall
    float t = curr_dist / (curr_dist - predict_dist);
    Vector2D intersection = pos + t * movement;

    // Check if the intersection point is within the wall segment bounds
    Vector2D to_intersection = intersection - this->start;
    float projection = Vector2D::dot(to_intersection, Vector2D::normalize(wall_vec));

    // Define epsilon width
    const float EPSILON = 1;

    // Ensure the intersection point lies within the wall segment with epsilon width
    return (projection >= -EPSILON && projection <= Vector2D::magnitude(wall_vec) + EPSILON);
}

void Wall::handleCollision(Particle& particle) const
{
    // Impulse
    Vector2D impulse = 2 * Vector2D::dot(-particle.getVelocity() * particle.getMass(), this->normal) * this->normal;

    // Apply impulse
    particle.applyImpulse(impulse);

}
