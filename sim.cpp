// sim.cpp
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <list>
#include <cmath>
#include <iostream>

#include "bondEnergies.h"

class Bond;

namespace Constants
{
    // Unit Conversions
    constexpr double PICO_TO_NANO = 0.001;

    // Simulation Parameters
    constexpr float TIMESTEP_PS = 0.005f;
    constexpr float SIMULATION_BOX_HALF_WIDTH = 1.0f;
    constexpr float DAMPING_FACTOR = 0.99995f; // Slight energy loss per step
    constexpr float WALL_COLLISION_DAMPING = -0.9f;

    // REFACTORED: Solver parameters
    constexpr int SOLVER_ITERATIONS = 5; 
    constexpr float CONSTRAINT_STRENGTH = 0.5f; 

    // Bonding Logic Parameters (Unchanged)
    constexpr float BOND_FORMATION_DISTANCE_SCALE = 1.1f;
    constexpr float BOND_BREAK_STRETCH_FACTOR = 1.8f;
    constexpr float INITIAL_SEPARATION_FACTOR = 1.1f;
}

enum class AtomType
{
    O, H, C, N, S, P, F, Cl, Br, I, He, Count
};

struct Atom
{
    float x, y;
    float old_x, old_y; 
    float ax, ay;       
    AtomType type;
    int valence;
    int covalentRadius;
    float atomicMass;
    std::vector<Bond*> bonds;
    void AddBond(Bond* bond) { bonds.push_back(bond); }
    void RemoveBond(Bond* bond) {
        auto it = std::remove(bonds.begin(), bonds.end(), bond);
        bonds.erase(it, bonds.end());
    }
};

class Bond
{
public:
    Atom* atom1;
    Atom* atom2;
    int order;
    double idealLength;
    double bondEnergy;

    Bond(Atom* a1, Atom* a2, int ord) : atom1(a1), atom2(a2), order(ord) {
        idealLength = static_cast<double>(a1->covalentRadius + a2->covalentRadius) * Constants::PICO_TO_NANO;
        const char* sym1 = BondData::atomTypeToSymbol(static_cast<int>(a1->type));
        const char* sym2 = BondData::atomTypeToSymbol(static_cast<int>(a2->type));
        bondEnergy = BondData::getBondEnthalpy(sym1, sym2, order);
    }

    float calculateLength() const {
        return std::sqrt(std::pow(atom2->x - atom1->x, 2) + std::pow(atom2->y - atom1->y, 2));
    }

    Atom* getOtherAtom(Atom* current) {
        return (current == atom1) ? atom2 : atom1;
    }
};

// --- Simulation Globals ---
extern "C"
{
const int MAX_ATOMS = 1000;
const int MAX_BONDS = MAX_ATOMS * 5;
std::array<Atom, static_cast<size_t>(AtomType::Count)> atomTemplates;
Atom atoms[MAX_ATOMS];
std::list<Bond> bondPool;
int atomCount = 0;
int simulation_step_counter = 0;

float bond_vertices[MAX_BONDS * 5]; 
int visible_bond_count = 0;

// --- Forward Declarations ---
void updatePositions(float dt);
void solveConstraints();
void applyBoundaryConditions();
void checkAndModifyBonds();
void updateBondVertexBuffer();
Bond* findBondBetween(Atom* a1, Atom* a2);
float getIdealAngleFor(Atom* centralAtom); 


float randf(float a, float b) {
    return a + (static_cast<float>(rand()) / RAND_MAX) * (b - a);
}

void init_atoms() {
    atomTemplates[(size_t)AtomType::O]  = {0,0,0,0,0,0, AtomType::O,  2, 66, 16.0f};
    atomTemplates[(size_t)AtomType::H]  = {0,0,0,0,0,0, AtomType::H,  1, 31, 1.008f};
    atomTemplates[(size_t)AtomType::C]  = {0,0,0,0,0,0, AtomType::C,  4, 73, 12.01f};
    atomTemplates[(size_t)AtomType::N]  = {0,0,0,0,0,0, AtomType::N,  3, 71, 14.01f};
    atomTemplates[(size_t)AtomType::S]  = {0,0,0,0,0,0, AtomType::S,  2, 105, 32.07f};
    atomTemplates[(size_t)AtomType::P]  = {0,0,0,0,0,0, AtomType::P,  3, 107, 30.97f};
    atomTemplates[(size_t)AtomType::F]  = {0,0,0,0,0,0, AtomType::F,  1, 57, 19.00f};
    atomTemplates[(size_t)AtomType::Cl] = {0,0,0,0,0,0, AtomType::Cl, 1, 102, 35.45f};
    atomTemplates[(size_t)AtomType::Br] = {0,0,0,0,0,0, AtomType::Br, 1, 120, 79.90f};
    atomTemplates[(size_t)AtomType::I]  = {0,0,0,0,0,0, AtomType::I,  1, 139, 126.90f};
    atomTemplates[(size_t)AtomType::He] = {0,0,0,0,0,0, AtomType::He, 0, 28, 4.0026f};

    atomCount = 50;
    bondPool.clear();
    
    for (int i = 0; i < atomCount; ++i) {
        int typeIndex = rand() % 2 + 1; 
        atoms[i] = atomTemplates[typeIndex];

        bool tooClose;
        do {
            tooClose = false;
            atoms[i].x = randf(-0.9f * Constants::SIMULATION_BOX_HALF_WIDTH, 0.9f * Constants::SIMULATION_BOX_HALF_WIDTH);
            atoms[i].y = randf(-0.9f * Constants::SIMULATION_BOX_HALF_WIDTH, 0.9f * Constants::SIMULATION_BOX_HALF_WIDTH);

            for (int j = 0; j < i; ++j) {
                float dx = atoms[i].x - atoms[j].x;
                float dy = atoms[i].y - atoms[j].y;
                float distSq = dx * dx + dy * dy;
                float min_dist = (atoms[i].covalentRadius + atoms[j].covalentRadius) * Constants::PICO_TO_NANO * Constants::INITIAL_SEPARATION_FACTOR;
                if (distSq < (min_dist * min_dist)) {
                    tooClose = true;
                    break;
                }
            }
        } while (tooClose);
        
        float vx = randf(-1.5f, 1.5f);
        float vy = randf(-1.5f, 1.5f);
        atoms[i].old_x = atoms[i].x - vx * Constants::TIMESTEP_PS;
        atoms[i].old_y = atoms[i].y - vy * Constants::TIMESTEP_PS;
        atoms[i].ax = 0.0f;
        atoms[i].ay = 0.0f;
    }
    simulation_step_counter = 0;
}


void spawnAtom(float x, float y, float intitalVelocityX, float intitalVelocityY, AtomType type) {
    if (atomCount >= MAX_ATOMS) return;

    Atom& newAtom = atoms[atomCount++];
    newAtom = atomTemplates[static_cast<size_t>(type)];
    newAtom.x = x;
    newAtom.y = y;
    newAtom.old_x = x - intitalVelocityX * Constants::TIMESTEP_PS;
    newAtom.old_y = y - intitalVelocityY * Constants::TIMESTEP_PS;
    newAtom.ax = 0.0f;
    newAtom.ay = 0.0f;
}



void simulate(float dt) {
    updatePositions(dt);
    solveConstraints();
    applyBoundaryConditions();

    if (simulation_step_counter % 20 == 0) {
        checkAndModifyBonds();
    }
    updateBondVertexBuffer(); 
    simulation_step_counter++;
}

void updatePositions(float dt) {
    for (int i = 0; i < atomCount; ++i) {
        // Velocity is implicitly (x - old_x)
        float vel_x = (atoms[i].x - atoms[i].old_x) * Constants::DAMPING_FACTOR;
        float vel_y = (atoms[i].y - atoms[i].old_y) * Constants::DAMPING_FACTOR;
        
        // Store current position
        atoms[i].old_x = atoms[i].x;
        atoms[i].old_y = atoms[i].y;
        
        // Perform Verlet integration
        atoms[i].x += vel_x + atoms[i].ax * dt * dt;
        atoms[i].y += vel_y + atoms[i].ay * dt * dt;

        // Reset acceleration
        atoms[i].ax = 0;
        atoms[i].ay = 0;
    }
}

void solveConstraints() {
    float delta_x[MAX_ATOMS] = {0};
    float delta_y[MAX_ATOMS] = {0};

    for (int k = 0; k < Constants::SOLVER_ITERATIONS; ++k) {
        for (const auto& bond : bondPool) {
            Atom* a1 = bond.atom1;
            Atom* a2 = bond.atom2;

            float dx = a2->x - a1->x;
            float dy = a2->y - a1->y;
            float dist = std::sqrt(dx * dx + dy * dy);
            if (dist < 1e-9f) { dist = 1e-9f; }

            float diff = (dist - static_cast<float>(bond.idealLength)) / dist;
            float response = diff * Constants::CONSTRAINT_STRENGTH;

            float totalMass = a1->atomicMass + a2->atomicMass;
            if (totalMass < 1e-9f) continue;
            float w1 = (a2->atomicMass / totalMass);
            float w2 = (a1->atomicMass / totalMass);

            a1->x += response * dx * w1;
            a1->y += response * dy * w1;
            a2->x -= response * dx * w2;
            a2->y -= response * dy * w2;
        }

        std::fill(delta_x, delta_x + atomCount, 0.0f);
        std::fill(delta_y, delta_y + atomCount, 0.0f);

        for (int i = 0; i < atomCount; ++i) {
            Atom* centralAtom = &atoms[i];
            if (centralAtom->bonds.size() < 2) continue;

            float idealAngleDeg = getIdealAngleFor(centralAtom);
            float idealAngleRad = idealAngleDeg * M_PI / 180.0f;

            for (size_t j = 0; j < centralAtom->bonds.size(); ++j) {
                for (size_t l = j + 1; l < centralAtom->bonds.size(); ++l) {
                    Bond* bond1 = centralAtom->bonds[j];
                    Bond* bond2 = centralAtom->bonds[l];

                    Atom* outerAtom1 = bond1->getOtherAtom(centralAtom);
                    Atom* outerAtom2 = bond2->getOtherAtom(centralAtom);

                    float L1 = bond1->idealLength;
                    float L2 = bond2->idealLength;
                    float idealDist = sqrtf(L1*L1 + L2*L2 - 2*L1*L2*cosf(idealAngleRad));

                    float dx = outerAtom2->x - outerAtom1->x;
                    float dy = outerAtom2->y - outerAtom1->y;
                    float dist = sqrtf(dx * dx + dy * dy);

                    if (dist > 1e-9f) {
                        float diff = (dist - idealDist) / dist;
                        float response = diff * Constants::CONSTRAINT_STRENGTH * 0.5f;

                        float totalMass = outerAtom1->atomicMass + outerAtom2->atomicMass;
                        if (totalMass < 1e-9f) continue;
                        float w1 = (outerAtom2->atomicMass / totalMass);
                        float w2 = (outerAtom1->atomicMass / totalMass);

                        delta_x[outerAtom1 - atoms] += response * dx * w1;
                        delta_y[outerAtom1 - atoms] += response * dy * w1;
                        delta_x[outerAtom2 - atoms] -= response * dx * w2;
                        delta_y[outerAtom2 - atoms] -= response * dy * w2;
                    }
                }
            }
        }

        for (int i = 0; i < atomCount; ++i) {
            atoms[i].x += delta_x[i];
            atoms[i].y += delta_y[i];
        }


        for (int i = 0; i < atomCount; ++i) {
            for (int j = i + 1; j < atomCount; ++j) {
                if (findBondBetween(&atoms[i], &atoms[j])) continue;

                float dx = atoms[j].x - atoms[i].x;
                float dy = atoms[j].y - atoms[i].y;
                float distSq = dx * dx + dy * dy;

                float collisionDist = (atoms[i].covalentRadius + atoms[j].covalentRadius) * Constants::PICO_TO_NANO;
                float collisionDistSq = collisionDist * collisionDist;

                if (distSq < collisionDistSq) {
                    float dist = std::sqrt(distSq);
                    if (dist < 1e-9f) { dist = 1e-9f; }

                    float diff = (dist - collisionDist) / dist;
                    float response = diff * Constants::CONSTRAINT_STRENGTH;

                    float totalMass = atoms[i].atomicMass + atoms[j].atomicMass;
                    float w1 = (atoms[j].atomicMass / totalMass);
                    float w2 = (atoms[i].atomicMass / totalMass);

                    atoms[i].x += response * dx * w1;
                    atoms[i].y += response * dy * w1;
                    atoms[j].x -= response * dx * w2;
                    atoms[j].y -= response * dy * w2;
                }
            }
        }
    }
}

void applyBoundaryConditions() {
    const float BOX_HW = Constants::SIMULATION_BOX_HALF_WIDTH;
    for (int i = 0; i < atomCount; ++i) {
        float vx = (atoms[i].x - atoms[i].old_x);
        float vy = (atoms[i].y - atoms[i].old_y);

        if (atoms[i].x < -BOX_HW) {
            atoms[i].x = -BOX_HW;
            atoms[i].old_x = atoms[i].x - vx * Constants::WALL_COLLISION_DAMPING;
        }
        if (atoms[i].x > BOX_HW) {
            atoms[i].x = BOX_HW;
            atoms[i].old_x = atoms[i].x - vx * Constants::WALL_COLLISION_DAMPING;
        }
        if (atoms[i].y < -BOX_HW) {
            atoms[i].y = -BOX_HW;
            atoms[i].old_y = atoms[i].y - vy * Constants::WALL_COLLISION_DAMPING;
        }
        if (atoms[i].y > BOX_HW) {
            atoms[i].y = BOX_HW;
            atoms[i].old_y = atoms[i].y - vy * Constants::WALL_COLLISION_DAMPING;
        }
    }
}


int calculateUsedValence(Atom* atom) {
    return std::accumulate(atom->bonds.begin(), atom->bonds.end(), 0,
        [](int sum, Bond* b) { return sum + b->order; });
}

Bond* findBondBetween(Atom* a1, Atom* a2) {
    for (auto* bond_ptr : a1->bonds) {
        if (bond_ptr->getOtherAtom(a1) == a2) {
            return bond_ptr;
        }
    }
    return nullptr;
}

void createBond(Atom* a1, Atom* a2, int order) {
    if (bondPool.size() >= MAX_BONDS) return;
    if (findBondBetween(a1, a2)) return;
    bondPool.emplace_back(a1, a2, order);
    Bond* newBond = &bondPool.back();
    if (newBond->bondEnergy > 1.0) { // Only form meaningful bonds
        a1->AddBond(newBond);
        a2->AddBond(newBond);
    } else {
        bondPool.pop_back();
    }
}

void breakBond(Bond* bond) {
    if (!bond || !bond->atom1) return;
    bond->atom1->RemoveBond(bond);
    bond->atom2->RemoveBond(bond);
    bond->atom1 = nullptr; // Mark for removal
    bond->atom2 = nullptr;
}

void checkAndModifyBonds() {
    // Bond formation / strengthening
    for (int i = 0; i < atomCount; ++i) {
        for (int j = i + 1; j < atomCount; ++j) {
            Atom* a1 = &atoms[i];
            Atom* a2 = &atoms[j];
            float dist = std::sqrt(std::pow(a2->x - a1->x, 2) + std::pow(a2->y - a1->y, 2));
            float idealDist = (a1->covalentRadius + a2->covalentRadius) * Constants::PICO_TO_NANO;

            if (dist < idealDist * Constants::BOND_FORMATION_DISTANCE_SCALE) {
                const char* sym1 = BondData::atomTypeToSymbol((int)a1->type);
                const char* sym2 = BondData::atomTypeToSymbol((int)a2->type);
                int availableValence1 = a1->valence - calculateUsedValence(a1);
                int availableValence2 = a2->valence - calculateUsedValence(a2);
                Bond* existingBond = findBondBetween(a1, a2);
                if (existingBond) {
                    if (existingBond->order < 3 && availableValence1 >= 1 && availableValence2 >= 1) {
                         double upgradedEnergy = BondData::getBondEnthalpy(sym1, sym2, existingBond->order + 1);
                         if (upgradedEnergy > existingBond->bondEnergy) {
                             existingBond->order++;
                             existingBond->bondEnergy = upgradedEnergy;
                         }
                    }
                } else {
                    if (availableValence1 >= 1 && availableValence2 >= 1) {
                        createBond(a1, a2, 1);
                    }
                }
            }
        }
    }

    std::vector<Bond*> bondsToBreak;
    for (auto& bond : bondPool) {
        if (bond.calculateLength() > bond.idealLength * Constants::BOND_BREAK_STRETCH_FACTOR) {
            bondsToBreak.push_back(&bond);
        }
    }
    for(auto* b : bondsToBreak) {
        breakBond(b);
    }
    
    bondPool.remove_if([](const Bond& b){ return b.atom1 == nullptr; });
}

void updateBondVertexBuffer() {
    visible_bond_count = 0;
    int i = 0;
    for (const auto& bond : bondPool) {
        if (visible_bond_count >= MAX_BONDS) break;
        bond_vertices[i++] = bond.atom1->x;
        bond_vertices[i++] = bond.atom1->y;
        bond_vertices[i++] = bond.atom2->x;
        bond_vertices[i++] = bond.atom2->y;
        bond_vertices[i++] = static_cast<float>(bond.order); 
        visible_bond_count++;
    }
}

float getIdealAngleFor(Atom* centralAtom) {
    int bondCount = centralAtom->bonds.size();
    if (bondCount <= 1) return 0.0f; // No angle to form
    if (bondCount == 2) return 180.0f; // Linear
    if (bondCount == 3) return 120.0f; // Trigonal Planar
    if (bondCount >= 4) return 109.5f; // Tetrahedral (projected onto 2D)
    return 0.0f;
}

// --- Getters for External Use ---
Atom* get_atoms() { return atoms; }
int get_atom_count() { return atomCount; }
int get_atom_size() { return sizeof(Atom); }
int get_atom_type_offset() { return offsetof(Atom, type); }
int get_atom_covalent_radius_offset() { return offsetof(Atom, covalentRadius); }

float* get_bond_vertices() { return bond_vertices; }
int get_visible_bond_count() { return visible_bond_count; }

} // extern "C"