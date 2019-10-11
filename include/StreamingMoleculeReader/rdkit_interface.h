//
// Created by Austin Clyde on 2019-10-06.
//

#ifndef RDKITSV_RDKIT_INTERFACE_H
#define RDKITSV_RDKIT_INTERFACE_H

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/Fingerprints/MorganFingerprints.h>
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
    ExplicitBitVect* getFingerPrint(std::string const &smi) {
        RDKit::ROMol *mol1 = nullptr;
        mol1 = RDKit::SmilesToMol(smi);

        if (mol1 == nullptr) {
            std::cerr << "ERROPROORORORO" << std::endl;
            return nullptr;
        }

        auto res = RDKit::MorganFingerprints::getFingerprintAsBitVect(*mol1, 3, 1024);
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
            auto res = RDKit::MorganFingerprints::getFingerprintAsBitVect(*mol1, 3, 1024);
            delete mol1;
            mol1 = nullptr;
            return std::make_pair(tmp, res);
        } else {
            return {};
        }
    }


    template<int dset_size>
    struct FastMinMax {
        std::vector<std::pair<ExplicitBitVect *, int>> pointers;
        int size;

        explicit FastMinMax(std::vector<std::string> const& smis_to_add)  {
            int i = 0;
            for(auto const& s : smis_to_add) {
                if (i >= dset_size) {
                    break;
                }
                auto fp = getFingerPrint(s);
                if (fp != nullptr) {
                    pointers.push_back(std::make_pair(fp, fp->getNumOnBits()));
                } else {
                    std::cerr << "one did noe work..." << std::endl;
                }
            }
            std::cout << "Loaded all the fprints" << std::endl;
            size = i;
        }

        inline float operator ()(ExplicitBitVect *mol)
        {
//            float min_sim = 0;
            float max_sim = 0;
            int bbits =  mol->getNumOnBits();
            for (int i = 0; i < size - 1; i++) {
                float ans = tanimotoSim(*(std::get<0>(pointers[i])), std::get<1>(pointers[i]), *mol, bbits);
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
