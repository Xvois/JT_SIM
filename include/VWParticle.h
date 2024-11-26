#ifndef VW_PARTICLE_H
#define VW_PARTICLE_H

#include "Vector2D.h"
#include "Particle.h"

class VWParticle : public Particle
{
private:
    float epsilon;       // Depth of the potential well
    float sigma;         // Characteristic distance for the potential

public:
    VWParticle() = default;
    VWParticle(Vector2D position, Vector2D velocity, float mass, float epsilon, float sigma);

    // Force computation for van der Waals interactions
    Vector2D computeForce(const Particle& other) const;

    // Interaction with other particles
    void interact(const Particle& particle);

    // Update particle position and velocity
    void update(float dt) override;

    using Particle::getPosition;
    using Particle::getVelocity;
    float getSigma() const { return sigma; }
    float getEpsilon() const { return epsilon; }
};

#endif // VW_PARTICLE_H