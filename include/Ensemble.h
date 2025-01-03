// Ensemble.h
#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <memory>

#include "Particle.h"
#include "Quad.h"
#include "Wall.h"

class Ensemble {
protected:
    std::vector<std::unique_ptr<Particle>> particles;
    Wall* bounds;
    int num_of_bounds;

public:
    Ensemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds);
    [[nodiscard]] static Vector2D collisionImpulse(Particle* p1, Particle* p2, float dt);
    virtual void iterateParticles(float dt);
    void addParticles(std::vector<std::unique_ptr<Particle>> newParticles);
    [[nodiscard]] bool isEmpty() const { return particles.empty(); }
    [[nodiscard]] std::vector<std::unique_ptr<Particle>> const& getParticles() const { return particles; }
    [[nodiscard]] double getTemperature() const;
    [[nodiscard]] double getTemperatureInRegion(const Quad& region) const;
    [[nodiscard]] double getPressureInRegion(const Quad& region) const;
    void cullNotInRegion(const Quad& region);
    void cullFastMovers(float maxSpeed);
};

#endif // ENSEMBLE_H