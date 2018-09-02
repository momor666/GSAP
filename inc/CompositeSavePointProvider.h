// Copyright (c) 2018 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
#ifndef PCOE_COMPOSITESAVEPOINTPROVIDER_h
#define PCOE_COMPOSITESAVEPOINTPROVIDER_h

#include "ISavePointProvider.h"

namespace PCOE {
    /**
     * A container for objects that implement ISavePointProvider that acts as a single 
     * Save Point Provider
     *
     * @author Christopher Teubert
     * @since 1.2
     **/
    class CompositeSavePointProvider : public ISavePointProvider {
    public:
        std::vector<ISavePointProvider *> providers;
        
        std::set<Message::time_point> getSavePts() override {
            std::set<Message::time_point> savePts;
            for (auto && provider : providers) {
                std::set<Message::time_point> tmp = provider->getSavePts();
                for (auto && elem : tmp) {
                    savePts.insert(elem);
                }
            }
            return savePts;
        }
        
        void add(ISavePointProvider * provider) {
            providers.push_back(provider);
        }
    };
}
#endif