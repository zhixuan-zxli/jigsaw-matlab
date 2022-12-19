/* 
jigsaw_mex.cpp

  mesh = jigsaw_mex(opts, geom, init, hfun);

  opts [in]
  geom [in]
  init [in][optional]
  hfun [in][optional]
  mesh [out]

See lib_jigsaw.h for details on the arguments. 
*/

#include <string>
#include "lib_jigsaw.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
protected:
  ArrayFactory factory;

protected:
  jigsaw_msh_t toMsh(const StructArray &matlabMesh) const {
    jigsaw_msh_t msh;
    jigsaw_init_msh_t(&msh);

    auto fields = matlabMesh.getFieldNames();
    std::vector<MATLABFieldIdentifier> fieldNames(fields.begin(), fields.end());
    CharArray mshID = matlabMesh[0]["mshID"];

    if(mshID.toAscii() == "EUCLIDEAN-MESH") {
      msh._flags = JIGSAW_EUCLIDEAN_MESH;

      // translate point array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("point"))
        != fieldNames.cend()) 
      {
        TypedArray<double> coord = matlabMesh[0]["point"][0]["coord"];
        auto coordDim = coord.getDimensions();
        // 2D point array
        if(coordDim[1] == 2) {
          jigsaw_alloc_vert2(&msh._vert2, coordDim[0]);
          for(size_t i = 0; i < coordDim[0]; ++i) {
            msh._vert2._data[i]._ppos[0] = coord[i][0];
            msh._vert2._data[i]._ppos[1] = coord[i][1];
          }
        // 3D point array
        } else if(coordDim[1] == 3) {
          jigsaw_alloc_vert3(&msh._vert3, coordDim[0]);
          for(size_t i = 0; i < coordDim[0]; ++i) {
            msh._vert3._data[i]._ppos[0] = coord[i][0];
            msh._vert3._data[i]._ppos[1] = coord[i][1];
            msh._vert3._data[i]._ppos[2] = coord[i][2];
          }
        }
      }

      // translate edge array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("edge2"))
        != fieldNames.cend()) 
      {
        TypedArray<double> index = matlabMesh[0]["edge2"][0]["index"];
        auto indexDim = index.getDimensions();
        jigsaw_alloc_edge2(&msh._edge2, indexDim[0]);
        for(size_t i = 0; i < indexDim[0]; ++i) {
          msh._edge2._data[i]._node[0] = index[i][0];
          msh._edge2._data[i]._node[1] = index[i][1];
        }
      }

      // translate triangle array      
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("tria3"))
        != fieldNames.cend()) 
      {
        TypedArray<double> index = matlabMesh[0]["tria3"][0]["index"];
        auto indexDim = index.getDimensions();
        jigsaw_alloc_tria3(&msh._tria3, indexDim[0]);
        for(size_t i = 0; i < indexDim[0]; ++i) {
          msh._tria3._data[i]._node[0] = index[i][0];
          msh._tria3._data[i]._node[1] = index[i][1];
          msh._tria3._data[i]._node[2] = index[i][2];
        }
      }

    } else if(mshID.toAscii() == "EUCLIDEAN-GRID") {
      msh._flags = JIGSAW_EUCLIDEAN_GRID;

      // translate grid array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("point"))
        != fieldNames.cend())
      {
        CellArray coord = matlabMesh[0]["point"][0]["coord"];
        auto nd = coord.getNumberOfElements();
        // x grid
        if(nd >= 1) {
          TypedArray<double> xgrid = coord[0];
          auto xgridDim = xgrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh._xgrid, xgridDim);
          for(int i = 0; i < xgridDim; ++i)
            msh._xgrid._data[i] = xgrid[i];
        }
        // y grid
        if(nd >= 2) {
          TypedArray<double> ygrid = coord[1];
          auto ygridDim = ygrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh._ygrid, ygridDim);
          for(int i = 0; i < ygridDim; ++i)
            msh._ygrid._data[i] = ygrid[i];
        }
        // z grid
        if(nd >= 3) {
          TypedArray<double> zgrid = coord[2];
          auto zgridDim = zgrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh._zgrid, zgridDim);
          for(int i = 0; i < zgridDim; ++i)
            msh._zgrid._data[i] = zgrid[i];
        }
      }

      // translate value array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("value"))
        != fieldNames.cend())
      {
        TypedArray<double> value = matlabMesh[0]["value"];
        auto valueDim = value.getNumberOfElements();
        jigsaw_alloc_flt32(&msh._value, valueDim);
        for(int i = 0; i < valueDim; ++i)
          msh._value._data[i] = value[i];
      }

    }

    return msh;
  }

public:
  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    checkArguments(outputs, inputs);

    auto geom = toMsh(inputs[1]);

    outputs[0] = factory.createScalar<double>(0.0);
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    if (inputs.size() < 2 || outputs.size() != 1) {
        matlabPtr->feval(u"error", 
            0, std::vector<matlab::data::Array>({ factory.createScalar("Invalid intput or output.") }));
    }
  }
};
