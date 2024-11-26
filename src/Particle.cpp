//
// Created by Sonny Parker on 21/11/2024.
//
#include "../include/Particle.h"

int Particle::next_id = 0;

void Particle::update(float const dt)
{
    this->position += velocity * dt;
}

void Particle::applyImpulse(const Vector2D& impulse)
{
    this->velocity += impulse / mass;
}

Vector2D Particle::predictPosition(const float dt) const
{
    return this->position + velocity * dt;
}