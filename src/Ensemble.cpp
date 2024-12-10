#include "../include/Ensemble.h"
#include "../include/VWParticle.h"
#include <iostream>
#include "../include/Constants.h"



Ensemble::Ensemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds)
    : particles(std::move(particles)), bounds(bounds), num_of_bounds(num_of_bounds) {}


Vector2D Ensemble::collisionImpulse(Particle* p1, Particle* p2, const float dt) {
    // Get current positions
    Vector2D p1_pos = p1->getPosition();
    Vector2D p2_pos = p2->getPosition();

    // Calculate relative velocity
    Vector2D v1 = p1->getVelocity();
    Vector2D v2 = p2->getVelocity();
    Vector2D relativeVelocity = v1 - v2;

    // Calculate normal vector
    Vector2D normal = Vector2D::normalize(p1_pos - p2_pos);

    // Calculate relative velocity in terms of the normal direction
    double velocityAlongNormal = Vector2D::dot(relativeVelocity, normal);

    // If velocities are separating, no impulse is needed
    if (velocityAlongNormal > 0) {
        return Vector2D{};
    }

    // Calculate impulse scalar for perfectly elastic collision (e = 1)
    double j = -(2 * velocityAlongNormal) / (1 / p1->getMass() + 1 / p2->getMass());

    // Apply impulse
    Vector2D impulse = j * normal;
    return impulse;
}


void Ensemble::iterateParticles(float dt)
{
    for (size_t i = 0; i < particles.size(); i++)
    {
        // Get a raw pointer to the current particle
        Particle* particle = particles[i].get();

        // Handle wall collisions for this particle
        for (int j = 0; j < num_of_bounds; j++)
        {
            if (bounds[j].checkCollision(*particle, dt))
            {
                bounds[j].handleCollision(*particle);
            }
        }

        // Update particle's position and velocity
        particle->update(dt);
    }
}

void Ensemble::addParticles(std::vector<std::unique_ptr<Particle>> newParticles)
{
    particles.insert(particles.end(), std::make_move_iterator(newParticles.begin()), std::make_move_iterator(newParticles.end()));
    newParticles.clear();
}


double Ensemble::getTemperature() const
{
    double sum = 0;
    for (const auto& particle : particles)
    {
        sum += particle->getKineticEnergy();
    }
    return sum / (particles.size() * K_b);
}

double Ensemble::getTemperatureInRegion(const Quad& region) const
{
    double sum = 0;
    int count = 0;
    for (const auto& particle : particles)
    {
        Vector2D pos = particle->getPosition();
        if (region.contains(pos.x, pos.y))
        {
            sum += particle->getKineticEnergy();
            count++;
        }
    }
    return sum / (count * K_b);
}

double Ensemble::getPressureInRegion(const Quad& region) const
{
    double V = region.getVolume() * SCALING * SCALING;
    double T = getTemperatureInRegion(region);
    double v = V * N_A / particles.size();
    int N = 0;
    for (const auto& particle : particles)
    {
        Vector2D pos = particle->getPosition();
        if (region.contains(pos.x, pos.y))
        {
            N++;
        }
    }
    const auto target = particles[0].get();
    // Check if particles is a VWParticle vector
    if (auto vwParticle = dynamic_cast<VWParticle*>(target))
    {
        // https://en.wikipedia.org/wiki/Van_der_Waals_equation
        const double sigma = vwParticle->getSigma();
        const double epsilon = vwParticle->getEpsilon();
        const double b = 4 * N_A * (4/3 * M_PI * pow(sigma/2, 3) );
        // reduced second virial coefficient??
        const int I = 1;
        const double a = I * N_A * epsilon * b;
        return (R * T)/(v - b) - a / pow(v, 2);
    }
    // Ideal gas law
    return N * K_b * T / V;
}

void Ensemble::cullNotInRegion(const Quad& region)
{
    std::erase_if(particles,
                  [&region](const std::unique_ptr<Particle>& particle) {
                      Vector2D pos = particle->getPosition();
                      return !region.contains(pos.x, pos.y);
                  });
}

void Ensemble::cullFastMovers(float maxSpeed)
{
    std::erase_if(particles,
                  [maxSpeed](const std::unique_ptr<Particle>& particle) {
                      return Vector2D::magnitude(particle->getVelocity()) > maxSpeed;
                  });
}
