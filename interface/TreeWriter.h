#ifndef ANALYSISTREEWRITER_H
#define ANALYSISTREEWRITER_H

#include <memory>
#include <map>
#include <string>
#include <iostream>

#include <TTree.h>
#include <TString.h>

class TreeWriter{
    public:
        TreeWriter();
        ~TreeWriter();

        void InitFloatBranch( const TString& name );
        void InitIntBranch(   const TString& name );

        void FillFloatBranch( const TString& name, float value );
        void FillIntBranch(   const TString& name, long  value );
        
        void InitFloatArray( const TString& name, const TString& lengthVariable );
        void InitIntArray(   const TString& name, const TString& lengthVariable );

        void FillFloatArray( const TString& name, int index, float value );
        void FillIntArray(   const TString& name, int index, long value );
    
        void SetTreeBranches(TTree* tree);

        void SetDefaultValues();
        void Dump() const;

    private:
        std::map<TString, Float_t> floatValues;
        std::map<TString, long>    intValues;

        float floatDefault = -9.;
        long  intDefault   = -1;

        std::map<TString, Float_t*> floatArrays;
        std::map<TString, long*>    intArrays;
        std::map<TString, TString>  arrayLength;

        int maxLength = 500;
        
};

#endif
