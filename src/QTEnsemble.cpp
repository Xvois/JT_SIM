//
// Created by Sonny Parker on 25/11/2024.
//
#include "../include/QTEnsemble.h"
#include "../include/VWParticle.h"
#include <iostream>
#include "../include/Constants.h"

QTEnsemble::QTEnsemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds, const Quad& boundary)
    : Ensemble(std::move(particles), bounds, num_of_bounds), tree(boundary) {
    for (auto& particle : this->particles) {
        tree.insert(particle.get());
    }
}

void QTEnsemble::iterateParticles(float dt) {
    // Clear and reinsert particles into the quadtree
    tree.clear();
    for (auto& particle : particles) {
        tree.insert(particle.get());
    }

    for (size_t i = 0; i < particles.size(); i++) {
        Particle* particle = particles[i].get();

        // Handle wall collisions for this particle
        for (int j = 0; j < num_of_bounds; j++) {
            if (bounds[j].checkCollision(*particle, dt)) {
                bounds[j].handleCollision(*particle);
            }
        }

        // Check for interactions with other particles
        if (auto vwParticle = dynamic_cast<VWParticle*>(particle)) {
            Quad range(particle->getPosition().x, particle->getPosition().y, 5, 5);
            std::vector<Particle*> neighbors;
            tree.query(range, neighbors);

            for (Particle* neighbor : neighbors) {
                if (neighbor != particle) {
                    vwParticle->interact(*neighbor);
                }
            }
        }

        // Update particle's position and velocity
        particle->update(dt);
    }
}

double QTEnsemble::getPressureInRegion(const Quad& region) const
{
    const double V = region.getVolume();
    int N = 0;
    for (const auto& particle : particles)
    {
        if (region.contains(particle->getPosition().x, particle->getPosition().y))
        {
            N++;
        }
    }

    const double v = V * N_A / N;
    const double T = getTemperatureInRegion(region);

    return 0;
}


void QTEnsemble::draw(sf::RenderWindow& window, bool showTree) const {
    Ensemble::draw(window);
    if (showTree)
    {
        tree.draw(window);
    }
}