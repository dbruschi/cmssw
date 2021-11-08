#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <vector>
#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>

template<typename T>
class FlattenedValueMapVectorTableProducer : public edm::stream::EDProducer<> {
    public:

        FlattenedValueMapVectorTableProducer( edm::ParameterSet const & params ):
            name_( params.getParameter<std::string>("name") ),
            doc_(params.existsAs<std::string>("doc") ? params.getParameter<std::string>("doc") : ""),
            src_(consumes<edm::View<T>>(params.getParameter<edm::InputTag>("src"))),
            cut_(params.existsAs<std::string>("cut") ? params.getParameter<std::string>("cut") : "", true)
        {
            edm::ParameterSet const & varsPSet = params.getParameter<edm::ParameterSet>("variables");
            for (const std::string & vname : varsPSet.getParameterNamesForType<edm::ParameterSet>()) {
                const auto & varPSet = varsPSet.getParameter<edm::ParameterSet>(vname);
                const std::string & type = varPSet.getParameter<std::string>("type");
                if (type == "std::vector<int>") {
                    intVecMaps_.emplace_back(consumes<edm::ValueMap<std::vector<int>>>(varPSet.getParameter<edm::InputTag>("src")));
                    intNames_.emplace_back(vname);
                    intDocs_.emplace_back(varPSet.getParameter<std::string>("doc"));
                    intPrecisions_.emplace_back(varPSet.getParameter<int>("precision"));
                }
                else if (type == "std::vector<float>") {
                    floatVecMaps_.emplace_back(consumes<edm::ValueMap<std::vector<float>>>(varPSet.getParameter<edm::InputTag>("src")));
                    floatNames_.emplace_back(vname);
                    floatDocs_.emplace_back(varPSet.getParameter<std::string>("doc"));
                    floatPrecisions_.emplace_back(varPSet.getParameter<int>("precision"));
                }
                else throw cms::Exception("Configuration", "unsupported type "+type+" for variable "+vname);
            }

            for (size_t i = 0 ; i < intVecMaps_.size(); i++) {
                produces<nanoaod::FlatTable>("intcounts"+std::to_string(i));
                produces<nanoaod::FlatTable>("intvec"+std::to_string(i));
            }
            for (size_t i = 0 ; i < floatVecMaps_.size(); i++) {
                produces<nanoaod::FlatTable>("floatcounts"+std::to_string(i));
                produces<nanoaod::FlatTable>("floatvec"+std::to_string(i));
            }
        }

        ~FlattenedValueMapVectorTableProducer() override {}

        template<typename P>
        std::vector<P> readVals(const edm::ValueMap<std::vector<P>>& vmap, edm::PtrVector<T>& objs, std::vector<int>& sizes) {
            std::vector<P> allvals;
            // Will be at least this long
            allvals.reserve(objs.size());
            for (size_t i = 0; i < objs.size(); i++) {
                auto& vals = vmap[objs[i]];
                sizes[i] = vals.size();
                allvals.insert(std::end(allvals), std::begin(vals), std::end(vals));
            }

            return allvals;
        }

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
            edm::Handle<edm::View<T>> src;
            iEvent.getByToken(src_, src);

            edm::PtrVector<T> objs;
            for (size_t i = 0; i < src->size(); ++i) {
                edm::Ptr<T> obj = src->ptrAt(i);
                if (cut_(*obj)) { 
                    objs.push_back(obj);
                }
            }

            std::vector<int> sizes(objs.size(), 0);
            for (size_t i = 0; i < intVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<int>>> vmap;
                iEvent.getByToken(intVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);

                auto intsizetab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_+intNames_[i]+"_Counts", false, false);
                intsizetab->template addColumn<int>("", sizes, "Number of entries per object", nanoaod::FlatTable::IntColumn, -1);
                intsizetab->setDoc(doc_);
                auto intvectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_+intNames_[i], false, false);
                intvectab->template addColumn<int>("Vals", results, intDocs_[i], nanoaod::FlatTable::IntColumn, intPrecisions_[i]);
                intvectab->setDoc(doc_);

                iEvent.put(std::move(intsizetab), "intcounts"+std::to_string(i));
                iEvent.put(std::move(intvectab), "intvec"+std::to_string(i));
            }
            std::fill(sizes.begin(), sizes.end(), 0);
            for (size_t i = 0; i < floatVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<float>>> vmap;
                iEvent.getByToken(floatVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);
                
                auto floatsizetab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_+floatNames_[i]+"_Counts", false, false);
                floatsizetab->template addColumn<int>("", sizes, "Number of entries per object", nanoaod::FlatTable::IntColumn, -1);
                floatsizetab->setDoc(doc_);
                auto floatvectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_+floatNames_[i], false, false);
                floatvectab->template addColumn<float>("Vals", results, floatDocs_[i], nanoaod::FlatTable::FloatColumn, floatPrecisions_[i]);
                floatvectab->setDoc(doc_);

                iEvent.put(std::move(floatsizetab), "floatcounts"+std::to_string(i));
                iEvent.put(std::move(floatvectab), "floatvec"+std::to_string(i));
            }
        }

    protected:
        const std::string name_; 
        const std::string doc_;
        const edm::EDGetTokenT<edm::View<T>> src_;
        const StringCutObjectSelector<T> cut_;
        std::vector<edm::EDGetTokenT<edm::ValueMap<std::vector<int>>>> intVecMaps_;
        std::vector<edm::EDGetTokenT<edm::ValueMap<std::vector<float>>>> floatVecMaps_;
        std::vector<std::string> intNames_;
        std::vector<std::string> floatNames_;
        std::vector<std::string> intDocs_;
        std::vector<std::string> floatDocs_;
        std::vector<int> intPrecisions_;
        std::vector<int> floatPrecisions_;
};
