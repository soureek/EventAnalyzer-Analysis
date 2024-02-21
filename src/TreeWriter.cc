#include "EventAnalyzer/Analysis/interface/TreeWriter.h"

using namespace std;

TreeWriter::TreeWriter() {}
TreeWriter::~TreeWriter() 
{
    for(std::map<TString, Float_t*>::iterator it=floatArrays.begin(); it!=floatArrays.end(); it++)
        delete[] it->second;

    for(std::map<TString, long*>::iterator it=intArrays.begin(); it!=intArrays.end(); it++)
        delete[] it->second;
}

void TreeWriter::InitFloatBranch( const TString& name )
{
    //floatValues[name] = 0.;
    floatValues[name] = (float)floatDefault;
}

void TreeWriter::InitIntBranch( const TString& name )
{
    //intValues[name] = 1.;
    intValues[name] = (int)intDefault;
}

void TreeWriter::InitFloatArray( const TString& name, const TString& lengthVariable )
{
    floatArrays[name] = new float[maxLength];
    for(int i = 0; i<maxLength; ++i)
    {
        floatArrays[name][i] = floatDefault;
    }
    arrayLength[name] = lengthVariable;
}

void TreeWriter::InitIntArray( const TString& name, const TString& lengthVariable )
{
    intArrays[name] = new long[maxLength];
    for(int i = 0; i<maxLength; ++i)
    {
        intArrays[name][i] = intDefault;
    }
    arrayLength[name] = lengthVariable;
}



void TreeWriter::FillFloatBranch( const TString& name, float value )
{
    floatValues[name] = value;
}

void TreeWriter::FillIntBranch( const TString& name, long value )
{
    intValues[name] = value;
}

void TreeWriter::FillFloatArray( const TString& name, int index, float value )
{
    floatArrays[name][index] = value;
}

void TreeWriter::FillIntArray( const TString& name, int index, long value )
{
    intArrays[name][index] = value;
}



void TreeWriter::SetTreeBranches(TTree* tree)
{
    auto iF = floatValues.begin();
    while(iF != floatValues.end())
    {
        tree->Branch(iF->first, &(iF->second), iF->first+"/F");
        iF++;
    }
    auto iI = intValues.begin();
    while(iI != intValues.end())
    {
        tree->Branch(iI->first, &(iI->second), iI->first+"/I");
        iI++;
    }
    auto iFA = floatArrays.begin();
    while(iFA != floatArrays.end())
    {
        tree->Branch(iFA->first, iFA->second, iFA->first+"["+arrayLength[iFA->first]+"]/F");
        iFA++;
    }
    auto iIA = intArrays.begin();
    while(iIA != intArrays.end())
    {
        tree->Branch(iIA->first, iIA->second, iIA->first+"["+arrayLength[iIA->first]+"]/I");
        iIA++;
    }

}

void TreeWriter::SetDefaultValues()
{
	
    auto iF = floatValues.begin();
    while(iF != floatValues.end())
    {
        iF->second = floatDefault;
        iF++;
    }
    auto iI = intValues.begin();
    while(iI != intValues.end())
    {
        iI->second = intDefault;
        iI++;
    }
    auto iFA = floatArrays.begin();
    while(iFA != floatArrays.end())
    {   
        for(int i = 0; i<maxLength; ++i)
        {
            iFA->second[i] = floatDefault;
        }
        iFA++;
    }
    auto iIA = intArrays.begin();
    while(iIA != intArrays.end())
    {
        for(int i = 0; i<maxLength; ++i)
        {
            iIA->second[i] = intDefault;
        }
        iIA++;
    }
}

void TreeWriter::Dump() const
{
    cout << "======== TEST EVENT ========" << endl;
    auto iF = floatValues.begin();
    while(iF != floatValues.end())
    {
        cout << iF->first << ": " << iF->second << endl;
        iF++;
    }
    auto iI = intValues.begin();
    while(iI != intValues.end())
    {
        cout << iI->first << ": " << iI->second << endl;
        iI++;
    }
    auto iFA = floatArrays.begin();
    while(iFA != floatArrays.end())
    {
        TString len = arrayLength.at(iFA->first);
        cout << iFA->first << ": ";
        for(int i = 0; i<intValues.at(len); ++i)
        {
            cout << iFA->second[i] << " ";
        }
        cout << endl;
        iFA++;
    }
    auto iIA = intArrays.begin();
    while(iIA != intArrays.end())
    {
        TString len = arrayLength.at(iIA->first);
        cout << iIA->first << ": ";
        for(int i = 0; i<intValues.at(len); ++i)
        {
            cout << iIA->second[i] << " ";
        }
        cout << endl;
        iIA++;
    }
    cout << "=============================" << endl;
}
