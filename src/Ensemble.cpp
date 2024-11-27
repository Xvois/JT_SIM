#include "../include/Ensemble.h"
#include "../include/VWParticle.h"
#include <iostream>
#include <SFML/Graphics/CircleShape.hpp>
#include "../include/Constants.h"



Ensemble::Ensemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds)
    : particles(std::move(particles)), bounds(bounds), num_of_bounds(num_of_bounds) {}

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

        // Handle wall collisions for this particle (again, to resolve corners)
        for (int j = 0; j < num_of_bounds; j++)
        {
            if (bounds[j].checkCollision(*particle, dt))
            {
                bounds[j].handleCollision(*particle);
            }
        }

        // Check for interactions with other particles
        if (auto vwParticle = dynamic_cast<VWParticle*>(particle))
        {
            for (size_t j = 0; j < particles.size(); j++)
            {
                if (i != j) // Avoid self-interaction
                {
                    vwParticle->interact(*particles[j]);
                }
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
    double V = region.getVolume();
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
        // UNSURE??
        const int I = 1;
        const double a = I * N_A * epsilon * b;
        return (R * T)/(v - b) - a / pow(v, 2);
    } else
    {
        // Ideal gas law
        return N * K_b * getTemperature() / V;
    }
}

void Ensemble::cullNotInRegion(const Quad& region)
{
    std::erase_if(particles,
                  [&region](const std::unique_ptr<Particle>& particle) {
                      Vector2D pos = particle->getPosition();
                      return !region.contains(pos.x, pos.y);
                  });
}

void Ensemble::draw(sf::RenderWindow& window) const
{
    for (const auto& particle : particles)
    {
        Vector2D pos = particle->getPosition();
        sf::CircleShape shape(0.5); // Radius of 1 pixel
        shape.setPosition(pos.x, pos.y);
        shape.setFillColor(sf::Color::White);
        window.draw(shape);
    }
    for (int i = 0; i < num_of_bounds; i++)
    {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(bounds[i].getStart().x, bounds[i].getStart().y)),
            sf::Vertex(sf::Vector2f(bounds[i].getEnd().x, bounds[i].getEnd().y))};
        window.draw(line, 2, sf::Lines);
    }
}