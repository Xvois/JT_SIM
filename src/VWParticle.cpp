//
// Created by Sonny Parker on 22/11/2024.
//
#include "../include/VWParticle.h"

#include <iostream>
#include <ostream>

VWParticle::VWParticle(Vector2D position, Vector2D velocity, float mass, float epsilon, float sigma)
    : Particle(position, velocity, mass), epsilon(epsilon), sigma(sigma) {}


// Compute the van der Waals force between this particle and another particle
Vector2D VWParticle::computeForce(const Particle& other) const {
    const Vector2D displacement = other.getPosition() - this->position;
    const double r = Vector2D::magnitude(displacement);

    const double LJ_force = - 24 * epsilon * (2 * pow(sigma, 12) / pow(r, 13) - pow(sigma, 6) / pow(r, 7));

    // Return the force as a vector
    return LJ_force * Vector2D::normalize(displacement);
}

void VWParticle::interact(const Particle& particle)
{
    acceleration += computeForce(particle) / mass;
}

// Update position and velocity based on net force and time step
void VWParticle::update(float dt) {
    // Verlet integration
    position += velocity * dt + 0.5f * acceleration * dt * dt;
    velocity += acceleration * dt;

    // ARGHHH!!! I FORGOT TO PUT THIS IN
    // AND HAVE SPENT EASILY >5HRS DEBUGGING IT!
    acceleration = Vector2D{};
}