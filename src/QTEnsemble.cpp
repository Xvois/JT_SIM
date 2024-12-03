//
// Created by Sonny Parker on 25/11/2024.
//
#include "../include/QTEnsemble.h"
#include "../include/VWParticle.h"
#include <iostream>
#include <thread>
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

    const unsigned int numThreads = std::thread::hardware_concurrency();
    const unsigned long chunkSize = particles.size() / numThreads;
    std::vector<std::thread> threads;

    auto updateChunk = [this, dt](int start, int end) {
        for (int i = start; i < end; ++i) {
            Particle* particle = particles[i].get();



            // Check for interactions with other particles
            if (auto vwParticle = dynamic_cast<VWParticle*>(particle)) {
                float sigma = vwParticle->getSigma();
                Quad range(vwParticle->getPosition().x, vwParticle->getPosition().y, 100 * sigma / SCALING, 100 * sigma / SCALING);
                std::vector<Particle*> neighbors;
                tree.query(range, neighbors);

                // Lennard-Jones interactions
                for (Particle* neighbor : neighbors) {
                    if (neighbor != particle) {
                        vwParticle->interact(*neighbor);
                    }
                }
            } else
            {
                // Handle particle collisions
                Quad range(particle->getPosition().x, particle->getPosition().y, 5, 5);
                std::vector<Particle*> neighbors;
                tree.query(range, neighbors);
                for (Particle* neighbor : neighbors) {
                    if (neighbor != particle) {
                        const Vector2D impulse = collisionImpulse(particle, neighbor, dt);
                        if ( (impulse.x != 0 || impulse.y != 0)) {
                            particle->applyImpulse(impulse);
                            neighbor->applyImpulse(-impulse);
                        }
                    }
                }
            }

            // Handle wall collisions for this particle
            for (int j = 0; j < num_of_bounds; j++) {
                if (bounds[j].checkCollision(*particle, dt)) {
                    bounds[j].handleCollision(*particle);
                }
            }

            // Update particle's position and velocity
            particle->update(dt);
        }
    };

    for (int i = 0; i < numThreads; ++i) {
        int start = i * chunkSize;
        int end = (i == numThreads - 1) ? particles.size() : start + chunkSize;
        threads.emplace_back(updateChunk, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }
}


void QTEnsemble::draw(sf::RenderWindow& window, bool showTree) const {
    Ensemble::draw(window);
    if (showTree)
    {
        tree.draw(window);
    }
}