//
// Created by Austin Clyde on 2019-10-06.
//

#ifndef RDKITSV_RDKIT_INTERFACE_H
#define RDKITSV_RDKIT_INTERFACE_H

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>

namespace SMR {

    boost::optional<std::string> getCannonicalSmileFromSmile(std::string const &smi) {
        RDKit::ROMol *mol1 = nullptr;
        try {
            mol1 = RDKit::SmilesToMol(smi);
        } catch (...) {
            if (mol1 != nullptr) {
                delete mol1;
                mol1 = nullptr;
            }
            return {};
        }

        if (mol1 != nullptr) {
            auto tmp = RDKit::MolToSmiles(*mol1);
            delete mol1;
            mol1 = nullptr;
            return tmp;
        } else {
            return {};
        }
    }
}

#endif //RDKITSV_RDKIT_INTERFACE_H
