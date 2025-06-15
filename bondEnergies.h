#pragma once
#include <unordered_map>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>


namespace BondData
{
    const char* atomTypeToSymbol(int type_enum) {
        static const std::vector<const char*> symbols = {"O", "H", "C", "N", "S", "P", "F", "Cl", "Br", "I", "He"};
        if (type_enum >= 0 && type_enum < symbols.size()) {
            return symbols[type_enum];
        }
        return "?";
    }


    const static std::unordered_map<std::string, double> BOND_DISSOCIATION_ENERGIES =
    {
        {"H-H", 436.0}, {"H-C", 337.2}, {"H-N", 314.0}, {"H-O", 428.0}, {"H-F", 568.6}, {"H-Cl", 431.8}, {"H-Br", 365.7}, {"H-I", 298.7},
        {"H-S", 344.0}, {"H-P", 343.0}, {"H-Si", 298.5},

        {"C-C", 368.0}, // H3C-CH3
        {"C=C", 682.0}, // H2C=CH2
        {"C#C", 962.0}, // HC#CH
        {"C-H", 431.0}, // in CH4
        {"C-O", 335.0}, // CH3-OCH3
        {"C=O", 732.0}, // H2C=O
        {"C#O", 1076.5},
        {"C-N", 331.0}, // CH3-NH2
        {"C#N", 937.0}, // HC#N
        {"C-Cl", 339.0},// CH3-Cl
        {"C-Br", 284.0},// CH3-Br
        {"C-I", 232.0},  // CH3-I
        {"C-F", 452.0},  // CH3-F
        {"C-S", 305.0},  // CH3-SH

        {"N-N", 167.0}, // in N2H4, approx value
        {"N=N", 456.0}, // HN=NH
        {"N#N", 945.3},
        {"N-O", 201.0}, // approx from various sources
        {"N-H", 391.0}, // in NH3
        {"N-Cl", 389.0},
        {"N-F", 301.0},

        {"O-O", 157.3}, // CH3O-OCH3
        {"O=O", 498.3},
        {"O-H", 498.7}, // H-OH
        {"O-F", 222.0},
        {"O-Cl", 272.0},
        {"O-Br", 235.1},
        {"O-I", 184.0},

        {"F-F", 156.9},
        {"Cl-Cl", 242.6},
        {"Br-Br", 193.9},
        {"I-I", 152.5},

        {"P-P", 490.0},
        {"P-H", 343.0},
        {"P-O", 596.6},
        {"P-S", 346.0},
        {"P-Cl", 289.0},

        {"S-S", 272.0}, // HS-SH
        {"S-H", 381.0}, // H-SH
        {"S-O", 521.7},
        {"S-Cl", 255.0},
    };

    double getBondEnthalpy(const std::string& atom1, const std::string& atom2, int order)
    {
        std::string key1, key2;
        std::string bond_symbol;

        switch (order)
        {
            case 1: bond_symbol = "-"; break;
            case 2: bond_symbol = "="; break;
            case 3: bond_symbol = "#"; break;
            default: return 0.0;
        }

        key1 = atom1 + bond_symbol + atom2;
        key2 = atom2 + bond_symbol + atom1;

        auto it = BOND_DISSOCIATION_ENERGIES.find(key1);
        if (it != BOND_DISSOCIATION_ENERGIES.end())
        {
            return it->second;
        }

        it = BOND_DISSOCIATION_ENERGIES.find(key2);
        if (it != BOND_DISSOCIATION_ENERGIES.end())
        {
            return it->second;
        }

        // Return 0 if no specific bond order is found.
        return 0.0;
    }
}