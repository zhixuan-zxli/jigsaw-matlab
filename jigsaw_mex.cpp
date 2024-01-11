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

#include <sstream>
#include <string>
#include "lib_jigsaw.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
protected:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
  ArrayFactory factory;
  std::ostringstream ostream;

  void displayOnMATLAB() {
    // Pass stream content to MATLAB fprintf function
    matlabPtr->feval(u"fprintf", 0,
        std::vector<Array>({ factory.createScalar(ostream.str()) }));
    // Clear stream buffer
    ostream.str("");
  }

  void displayDim(const ArrayDimensions &dim) {
    ostream << "[";
    for(auto n : dim)
      ostream << n << ",";
    ostream << "]";
  }

protected:
  jigsaw_msh_t * to_msh(const StructArray &matlabMesh) {
    jigsaw_msh_t * msh = new jigsaw_msh_t;
    jigsaw_init_msh_t(msh);

    auto fields = matlabMesh.getFieldNames();
    std::vector<MATLABFieldIdentifier> fieldNames(fields.begin(), fields.end());
    const CharArray mshID = matlabMesh[0]["mshID"];

    if(mshID.toAscii() == "EUCLIDEAN-MESH") {
      msh->_flags = JIGSAW_EUCLIDEAN_MESH;

      // translate point array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("point"))
        != fieldNames.cend()) 
      {
        const StructArray point = matlabMesh[0]["point"];
        const TypedArray<double> coord = point[0]["coord"];
        auto coordDim = coord.getDimensions();
        // 2D point array
        if(coordDim[1] == 3) {
          jigsaw_alloc_vert2(&msh->_vert2, coordDim[0]);
          for(size_t i = 0; i < coordDim[0]; ++i) {
            msh->_vert2._data[i]._ppos[0] = coord[i][0];
            msh->_vert2._data[i]._ppos[1] = coord[i][1];
            msh->_vert2._data[i]._itag = coord[i][2];
          }
        // 3D point array
        } else if(coordDim[1] == 4) {
          jigsaw_alloc_vert3(&msh->_vert3, coordDim[0]);
          for(size_t i = 0; i < coordDim[0]; ++i) {
            msh->_vert3._data[i]._ppos[0] = coord[i][0];
            msh->_vert3._data[i]._ppos[1] = coord[i][1];
            msh->_vert3._data[i]._ppos[2] = coord[i][2];
            msh->_vert3._data[i]._itag = coord[i][3];
          }
        }
      }

      // translate edge array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("edge2"))
        != fieldNames.cend()) 
      {
        const StructArray edge2 = matlabMesh[0]["edge2"];
        const TypedArray<double> index = edge2[0]["index"];
        auto indexDim = index.getDimensions();
        jigsaw_alloc_edge2(&msh->_edge2, indexDim[0]);
        for(size_t i = 0; i < indexDim[0]; ++i) {
          msh->_edge2._data[i]._node[0] = index[i][0] - 1;
          msh->_edge2._data[i]._node[1] = index[i][1] - 1;
          msh->_edge2._data[i]._itag = index[i][2];
        }
      }

      // translate bound array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("bound"))
        != fieldNames.cend()) 
      {
        const StructArray bound = matlabMesh[0]["bound"];
        const TypedArray<double> index = bound[0]["index"];
        auto indexDim = index.getDimensions();
        jigsaw_alloc_bound(&msh->_bound, indexDim[0]);
        for(size_t i = 0; i < indexDim[0]; ++i) {
          msh->_bound._data[i]._itag = index[i][0];
          msh->_bound._data[i]._indx = index[i][1] - 1;
          msh->_bound._data[i]._kind = index[i][2];
        }
      }

      // translate triangle array      
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("tria3"))
        != fieldNames.cend()) 
      {
        const StructArray tria3 = matlabMesh[0]["tria3"];
        const TypedArray<double> index = tria3[0]["index"];
        auto indexDim = index.getDimensions();
        jigsaw_alloc_tria3(&msh->_tria3, indexDim[0]);
        for(size_t i = 0; i < indexDim[0]; ++i) {
          msh->_tria3._data[i]._node[0] = index[i][0] - 1;
          msh->_tria3._data[i]._node[1] = index[i][1] - 1;
          msh->_tria3._data[i]._node[2] = index[i][2] - 1;
          msh->_tria3._data[i]._itag = index[i][3];
        }
      }

    } else if(mshID.toAscii() == "EUCLIDEAN-GRID") {
      msh->_flags = JIGSAW_EUCLIDEAN_GRID;

      // translate grid array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("point"))
        != fieldNames.cend())
      {
        const StructArray point = matlabMesh[0]["point"];
        const CellArray coord = point[0]["coord"];
        auto nd = coord.getNumberOfElements();
        // x grid
        if(nd >= 1) {
          const TypedArray<double> xgrid = coord[0];
          auto xgridDim = xgrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh->_xgrid, xgridDim);
          for(int i = 0; i < xgridDim; ++i)
            msh->_xgrid._data[i] = xgrid[i];
        }
        // y grid
        if(nd >= 2) {
          const TypedArray<double> ygrid = coord[1];
          auto ygridDim = ygrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh->_ygrid, ygridDim);
          for(int i = 0; i < ygridDim; ++i)
            msh->_ygrid._data[i] = ygrid[i];
        }
        // z grid
        if(nd >= 3) {
          const TypedArray<double> zgrid = coord[2];
          auto zgridDim = zgrid.getNumberOfElements();
          jigsaw_alloc_reals(&msh->_zgrid, zgridDim);
          for(int i = 0; i < zgridDim; ++i)
            msh->_zgrid._data[i] = zgrid[i];
        }
      }

      // translate value array
      if(std::find(fieldNames.cbegin(), fieldNames.cend(), MATLABFieldIdentifier("value"))
        != fieldNames.cend())
      {
        const TypedArray<double> value = matlabMesh[0]["value"];
        auto valueDim = value.getNumberOfElements();
        jigsaw_alloc_flt32(&msh->_value, valueDim);
        for(int i = 0; i < valueDim; ++i)
          msh->_value._data[i] = value[i];
      }

    }

    return msh;
  }

  StructArray from_msh(const jigsaw_msh_t *msh) {
    std::vector<std::string> fieldNames { "mshID" };
    if(msh->_vert2._data || msh->_vert3._data)
      fieldNames.push_back("point");
    if(msh->_edge2._data)
      fieldNames.push_back("edge2");
    if(msh->_tria3._data)
      fieldNames.push_back("tria3");

    StructArray matlabMesh = factory.createStructArray({1}, fieldNames);
    // process mshID
    if(msh->_flags == JIGSAW_EUCLIDEAN_MESH)
      matlabMesh[0]["mshID"] = factory.createCharArray("EUCLIDEAN-MESH");

    // process points
    if(msh->_vert2._data) {
      matlabMesh[0]["point"] = factory.createStructArray({1}, {"coord"});
      TypedArrayRef<Struct> point = matlabMesh[0]["point"];
      ArrayDimensions dim { msh->_vert2._size, 3 };
      point[0]["coord"] = factory.createArray<double>(dim);
      TypedArrayRef<double> coord = point[0]["coord"];
      for(size_t i = 0; i < dim[0]; ++i) {
        coord[i][0] = msh->_vert2._data[i]._ppos[0];
        coord[i][1] = msh->_vert2._data[i]._ppos[1];
        coord[i][2] = msh->_vert2._data[i]._itag;
      }
    } else if(msh->_vert3._data) {
      matlabMesh[0]["point"] = factory.createStructArray({1}, {"coord"});
      TypedArrayRef<Struct> point = matlabMesh[0]["point"];
      ArrayDimensions dim { msh->_vert3._size, 4 };
      point[0]["coord"] = factory.createArray<double>(dim);
      TypedArrayRef<double> coord = point[0]["coord"];
      for(size_t i = 0; i < dim[0]; ++i) {
        coord[i][0] = msh->_vert3._data[i]._ppos[0];
        coord[i][1] = msh->_vert3._data[i]._ppos[1];
        coord[i][2] = msh->_vert3._data[i]._ppos[2];
        coord[i][3] = msh->_vert3._data[i]._itag;
      }
    }

    // process edges
    if(msh->_edge2._data) {
      matlabMesh[0]["edge2"] = factory.createStructArray({1}, {"index"});
      TypedArrayRef<Struct> edge2 = matlabMesh[0]["edge2"];
      ArrayDimensions dim { msh->_edge2._size, 3 };
      edge2[0]["index"] = factory.createArray<double>(dim);
      TypedArrayRef<double> index = edge2[0]["index"];
      for(size_t i = 0; i < dim[0]; ++i) {
        index[i][0] = msh->_edge2._data[i]._node[0] + 1;
        index[i][1] = msh->_edge2._data[i]._node[1] + 1;
        index[i][2] = msh->_edge2._data[i]._itag;
      }
    }

    // process triangles
    if(msh->_tria3._data) {
      matlabMesh[0]["tria3"] = factory.createStructArray({1}, {"index"});
      TypedArrayRef<Struct> tria3 = matlabMesh[0]["tria3"];
      ArrayDimensions dim { msh->_tria3._size, 4 };
      tria3[0]["index"] = factory.createArray<double>(dim);
      TypedArrayRef<double> index = tria3[0]["index"];
      for(size_t i = 0; i < dim[0]; ++i) {
        index[i][0] = msh->_tria3._data[i]._node[0] + 1;
        index[i][1] = msh->_tria3._data[i]._node[1] + 1;
        index[i][2] = msh->_tria3._data[i]._node[2] + 1;
        index[i][3] = msh->_tria3._data[i]._itag;
      }
    }

    return matlabMesh;
  }

  jigsaw_jig_t to_jig(const StructArray &matlabJig)
  {
    jigsaw_jig_t jig;
    jigsaw_init_jig_t(&jig);

    auto fields = matlabJig.getFieldNames();
    std::vector<MATLABFieldIdentifier> fieldNames(fields.begin(), fields.end());

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("mesh_kern"))
      != fieldNames.end()) {
      const CharArray mesh_kern = matlabJig[0]["mesh_kern"];
      if(mesh_kern.toAscii() == "delfront")
        jig._mesh_kern = JIGSAW_KERN_DELFRONT;
      else if(mesh_kern.toAscii() == "delaunay")
        jig._mesh_kern = JIGSAW_KERN_DELAUNAY;
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("mesh_dims"))
      != fieldNames.end()) {
      const TypedArray<double> mesh_dims = matlabJig[0]["mesh_dims"];
      jig._mesh_dims = mesh_dims[0];
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("optm_qlim"))
      != fieldNames.end()) {
      const TypedArray<double> optm_qlim = matlabJig[0]["optm_qlim"];
      jig._optm_qlim = optm_qlim[0];
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("hfun_scal"))
      != fieldNames.end()) {
      const CharArray hfun_scal = matlabJig[0]["hfun_scal"];
      if(hfun_scal.toAscii() == "relative")
        jig._hfun_scal = JIGSAW_HFUN_RELATIVE;
      else if(hfun_scal.toAscii() == "absolute")
        jig._hfun_scal = JIGSAW_HFUN_ABSOLUTE;
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("hfun_hmax"))
      != fieldNames.end()) {
      const TypedArray<double> hfun_hmax = matlabJig[0]["hfun_hmax"];
      jig._hfun_hmax = hfun_hmax[0];
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("hfun_hmin"))
      != fieldNames.end()) {
      const TypedArray<double> hfun_hmin = matlabJig[0]["hfun_hmin"];
      jig._hfun_hmin = hfun_hmin[0];
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("geom_feat"))
      != fieldNames.end()) {
      const TypedArray<bool> geom_feat = matlabJig[0]["geom_feat"];
      jig._geom_feat = geom_feat[0];
    }

    if(std::find(fieldNames.begin(), fieldNames.end(), MATLABFieldIdentifier("mesh_top1"))
      != fieldNames.end()) {
      const TypedArray<bool> mesh_top1 = matlabJig[0]["mesh_top1"];
      jig._mesh_top1 = mesh_top1[0];
    }

    return jig;
  }

public:
  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    checkArguments(outputs, inputs);

    auto jig = to_jig(inputs[0]);

    jigsaw_msh_t * geom = to_msh(inputs[1]);
    jigsaw_msh_t * init = (inputs.size() >= 3 && inputs[2].getNumberOfElements() >= 1) ? to_msh(inputs[2]) : nullptr;
    jigsaw_msh_t * hfun = (inputs.size() >= 4 && inputs[3].getNumberOfElements() >= 1) ? to_msh(inputs[3]) : nullptr;

    // outputs[0] = from_msh(geom);
    jigsaw_msh_t mesh;
    jigsaw_init_msh_t(&mesh);
    auto ret = jigsaw(&jig, geom, init, hfun, &mesh);

    if(outputs.size() >= 1)
      outputs[0] = from_msh(&mesh);
    if(outputs.size() >= 2)
      outputs[1] = factory.createScalar<double>(ret);

    jigsaw_free_msh_t(&mesh);
    jigsaw_free_msh_t(geom); delete geom;
    if(init) { jigsaw_free_msh_t(init); delete init; }
    if(hfun) { jigsaw_free_msh_t(hfun); delete hfun; }
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    if (inputs.size() < 2 || outputs.size() > 2) {
        matlabPtr->feval(u"error", 
            0, std::vector<matlab::data::Array>({ factory.createScalar("Invalid intput or output.") }));
    }
  }
};
