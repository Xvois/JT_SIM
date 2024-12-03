//
// Created by Sonny Parker on 22/11/2024.
//
#include "../include/VWParticle.h"

#include <iostream>
#include <ostream>

#include "../include/Constants.h"

VWParticle::VWParticle(Vector2D position, Vector2D velocity, float mass, float epsilon, float sigma)
    : Particle(position, velocity, mass), epsilon(epsilon), sigma(sigma) {}


// Compute the van der Waals force between this particle and another particle
Vector2D VWParticle::computeForce(const Particle& other) const {
    const Vector2D displacement = other.getPosition() - this->position;
    const double r = Vector2D::magnitude(displacement) * SCALING + sigma/4;
    double LJ_force = -24 * epsilon * (2 * pow(sigma, 12) / pow(r, 13) - pow(sigma, 6) / pow(r, 7));

    // REMOVE SINGULARITY
    if (LJ_force > 2e-21) {
        LJ_force = 2e-21;
    }

    return LJ_force * Vector2D::normalize(displacement);
}

void VWParticle::interact(const Particle& particle)
{
    acceleration += computeForce(particle) / mass;
}

// Update position and velocity based on net force and time step
void VWParticle::update(float dt) {
    // Verlet integration
    position += (velocity / SCALING) * dt + 0.5f * (acceleration / SCALING) * dt * dt;
    velocity += (acceleration / SCALING) * dt;

    // ARGHHH!!! I FORGOT TO PUT THIS IN
    // AND HAVE SPENT EASILY >5HRS DEBUGGING IT!
    acceleration = Vector2D{};
}