//
// Created by Sonny Parker on 25/11/2024.
//

#ifndef QTENSEMBLE_H
#define QTENSEMBLE_H

#include <vector>
#include "Particle.h"
#include "Wall.h"
#include "QuadTree.h"
#include "Ensemble.h"

class QTEnsemble : public Ensemble {
private:
    QuadTree tree;

public:
    QTEnsemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds, const Quad& boundary);
    void iterateParticles(float dt) override;

    [[nodiscard]] QuadTree& getTree() { return tree; }

    using Ensemble::getPressureInRegion;
    using Ensemble::getTemperatureInRegion;
    using Ensemble::getTemperature;

};

#endif // QTENSEMBLE_H
