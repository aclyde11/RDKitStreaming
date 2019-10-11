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
#include <RDGeneral/export.h>
#include <DataStructs/ExplicitBitVect.h>

namespace SMR {

    template<typename X>
    inline float tanimotoSim(X const &A, X const &B) {
        int both = (A & B).getNumOnBits();
        int denom = A.getNumOnBits() + B.getNumOnBits() - both;
        return static_cast<float>(both) / static_cast<float>(denom);
    }

    template<typename X>
    inline float tanimotoSim(X const &A, int ac, X const &B) {
        int both = (A & B).getNumOnBits();
        int denom = A.getNumOnBits() + ac - both;
        return static_cast<float>(both) / static_cast<float>(denom);
    }

    template<typename X>
    inline float tanimotoSim(X const &A, int ac, X const &B, int bc) {
        int both = (A & B).getNumOnBits();
        int denom = bc + ac - both;
        return static_cast<float>(both) / static_cast<float>(denom);
    }


    //MUST BE VALID.
    ExplicitBitVect *getFingerPrint(std::string const &smi) {
        RDKit::ROMol *mol1 = nullptr;
        mol1 = RDKit::SmilesToMol(smi);
        auto res = RDKit::RDKFingerprintMol(*mol1, 1, 7, 1024);
        delete mol1;
        mol1 = nullptr;
        return res;
    }

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

    std::pair<boost::optional<std::string>, ExplicitBitVect *> getCannonicalSmileFromSmileFP(std::string const &smi) {
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
            auto res = RDKit::RDKFingerprintMol(*mol1, 1, 7, 1024);
            delete mol1;
            mol1 = nullptr;
            return std::make_pair(tmp, res);
        } else {
            return {};
        }
    }


    template<int dset_size>
    struct FastMinMax {
        ExplicitBitVect  *pointers[dset_size];
        int count[dset_size];

        FastMinMax(std::vector<std::string> const& smis_to_add) : pointers{nullptr} , count{0} {
            int i = 0;
            for(auto const& s : smis_to_add) {
                if (i >= dset_size) {
                    break;
                }
                auto fp = getFingerPrint(s);
                count[i] = fp->getNumOnBits();
                pointers[i] = fp;
                i++;
            }
        }

        inline float operator ()(ExplicitBitVect *mol)
        {
//            float min_sim = 0;
            float max_sim = 0;
            int bbits =  mol->getNumOnBits();
            for (int i = 0; i < dset_size; i++) {
                float ans = tanimotoSim(*(pointers[i]), count[i], *mol, bbits);
                if (ans > max_sim ) {
                    max_sim = ans;
                }
            }

            delete mol;
            mol = nullptr;
            return max_sim;
        }
    };
}

#endif //RDKITSV_RDKIT_INTERFACE_H
