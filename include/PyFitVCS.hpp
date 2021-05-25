#ifndef _PY_FIT_VCS_HPP_
#define _PY_FIT_VCS_HPP_
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <functional>
#include <vector>
#include <string>
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include "Integration.hpp"
#include "KuraevFadin.hpp"
#include "PyUtils.hpp"

typedef struct {
  PyObject_HEAD
  bool energySpread;
  double threshold;
  PyArrayObject* energy;
  PyArrayObject* energyErr;
  PyArrayObject* vcs;
  PyArrayObject* vcsErr;
  PyObject* bcsModelFCN;
  PyObject* effFCN;
  std::function<double(double, double)> effLambda;
  bool energy_spread;
  double errordef;
  PyObject* param_names;
} PyFitVCSObject;

static PyMemberDef PyFitVCS_members[] = {
  {"errordef", T_DOUBLE, offset_of(PyFitVCSObject, errordef), 0,
   "errordef"},
  {NULL}
};

static void
PyFitVCS_dealloc(PyFitVCSObject *self)
{
  Py_XDECREF(self->energy);
  Py_XDECREF(self->energyErr);
  Py_XDECREF(self->vcs);
  Py_XDECREF(self->vcsErr);
  Py_XDECREF(self->bcsModelFCN);
  Py_XDECREF(self->effFCN);
  Py_XDECREF(self->param_names);
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
PyFitVCS_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyFitVCSObject *self;
  self = reinterpret_cast<PyFitVCSObject*>(type->tp_alloc(type, 0));
  if (self != NULL) {
    self->energySpread = false;
    self->energy = NULL;
    self->energyErr = NULL;
    self->vcs = NULL;
    self->vcsErr = NULL;
    self->bcsModelFCN = NULL;
    self->effFCN = NULL;
    self->effLambda =
        [self](double x, double en) {
          PyObject *arglist = Py_BuildValue("(dd)", x, en);
          PyObject *rv = PyObject_CallObject(self->effFCN, arglist);
          double result = PyFloat_AS_DOUBLE(rv);
          Py_CLEAR(rv);
          Py_CLEAR(arglist);
          return result;
        };
  }
  return reinterpret_cast<PyObject*>(self);
}

static int
PyFitVCS_init(PyFitVCSObject *self, PyObject *args, PyObject *kwds) {
  static const char *kwlist[] =
      {"threshold", "energy", "vcs", "energy_err",
       "vcs_err", "bcs_model", "efficiency", NULL};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "dO!O!O!O!OO",
          const_cast<char**>(kwlist),
          &(self->threshold),
          &PyArray_Type, &(self->energy),
          &PyArray_Type, &(self->vcs),
          &PyArray_Type, &(self->energyErr),
          &PyArray_Type, &(self->vcsErr),
          &(self->bcsModelFCN), &(self->effFCN))) {
    return -1;
  }
  if (!PyCallable_Check(self->bcsModelFCN)) {
    PyErr_SetString(PyExc_TypeError, "IterISRSolverUseVCSFit: a callable Born cross section object is required");
    return -1;
  }
  if (!PyCallable_Check(self->effFCN)) {
    PyErr_SetString(PyExc_TypeError, "IterISRSolverUseVCSFit: a callable efficiency object is required");
    return -1;
  }
  PyObject* const inspect_module_name = PyUnicode_DecodeFSDefault("inspect");
  PyObject* const inspect_module = PyImport_Import(inspect_module_name);
  PyObject* const getargspec_function = PyObject_GetAttrString(inspect_module, "getfullargspec");
  PyObject* const argspec_call_args = PyTuple_New(1);
  PyTuple_SetItem(argspec_call_args, 0, self->bcsModelFCN);
  PyObject* const argspec = PyObject_CallObject(getargspec_function, argspec_call_args);
  //self->param_names = PyObject_GetAttrString(argspec, "args");
  // Py_XINCREF(self->param_names);
  PyObject* const f_args = PyObject_GetAttrString(argspec, "args");
  Py_ssize_t const num_args = PyList_Size(f_args);
  self->param_names = PyTuple_New(num_args - 1);
  // self->par_names.reserve(num_args - 1);
  for (Py_ssize_t i = 1; i < num_args; ++i)
  {
    PyObject* const arg = PyList_GetItem(f_args, i);
    Py_XINCREF(arg);
    PyTuple_SET_ITEM(self->param_names, i - 1, arg);
  }
  self->energy_spread = false;
  self->errordef = 1.;

  Py_XDECREF(f_args);
  Py_XINCREF(self->energy);
  Py_XINCREF(self->vcs);
  Py_XINCREF(self->energyErr);
  Py_XINCREF(self->vcsErr);
  Py_XINCREF(self->bcsModelFCN);
  Py_XINCREF(self->effFCN);
  return 0;
}

static PyObject *PyFitVCS_call(PyObject *callable, PyObject *args, PyObject* kwds) {
  auto self = reinterpret_cast<PyFitVCSObject*>(callable);
  Py_XINCREF(args);
  double* energyC = (double*) PyArray_DATA(self->energy);
  double* energyErrC = (double*) PyArray_DATA(self->energyErr);
  double* vcsC = (double*) PyArray_DATA(self->vcs);
  double* vcsErrC = (double*) PyArray_DATA(self->vcsErr);
  // !!! dimension check
  npy_intp *dims = PyArray_DIMS(self->energy);
  npy_intp dim = dims[0];
  double result = 0;
  Py_ssize_t paramsSize = PyTuple_Size(args);
  PyObject* argtuple = PyTuple_New(paramsSize + 1);
  PyTuple_SET_ITEM(argtuple, 0, PyFloat_FromDouble(0.));
  for (Py_ssize_t i = 0; i < paramsSize; ++i) {
    PyObject* obj = PyTuple_GET_ITEM(args, i);
    Py_XINCREF(obj);
    PyTuple_SET_ITEM(argtuple, i + 1, obj);
  }
  std::function<double(double)> bcsModelLambda =
      [&argtuple, self, kwds](double en) {
        PyTuple_SET_ITEM(argtuple, 0, PyFloat_FromDouble(en));
        PyObject* rv = PyObject_Call(self->bcsModelFCN, argtuple, kwds);
        double result = PyFloat_AS_DOUBLE(rv);
        Py_XDECREF(rv);
        return result;
      };
  for (npy_intp i = 0; i < dim; ++i) {
    double sT = self->threshold * self->threshold;
    const double energy = energyC[i];
    double modelVCS = 0;
    if (self->energy_spread) {
      std::function<double(double)> fcnVCSNoSpread =
          [bcsModelLambda, sT, self](double en) {
            const double s = en * en;
            double result = 0;
            if (s <= sT) {
              return result;
            }
            result = convolutionKuraevFadin(
                en, bcsModelLambda,
                0, 1. - sT / s,
                self->effLambda);
            return result;
          };
      const double sigmaEn2 = energyErrC[i] * energyErrC[i];
      modelVCS = gaussian_conv(energy, sigmaEn2, fcnVCSNoSpread);
    } else {
      modelVCS = convolutionKuraevFadin(
          energy, bcsModelLambda,
          0, 1. - sT / energy / energy,
          self->effLambda);
    }
    double dchi2 = (vcsC[i] - modelVCS) / vcsErrC[i];
    dchi2 *= dchi2;
    result += dchi2;
  }
  Py_XDECREF(argtuple);
  Py_XDECREF(args);
  return PyFloat_FromDouble(result);
}

static PyObject *PyFitVCS_minuit(PyFitVCSObject *self, PyObject *args, PyObject *kwds) {
  PyObject* const minuit_module_name = PyUnicode_DecodeFSDefault("iminuit");
  PyObject* const minuit_module = PyImport_Import(minuit_module_name);
  PyObject* self_obj = reinterpret_cast<PyObject*>(self);
  PyObject* const minuit_function = PyObject_GetAttrString(minuit_module, "Minuit");
  Py_ssize_t paramsSize = PyTuple_Size(args);
  PyObject* argtuple = PyTuple_New(paramsSize + 1);
  PyTuple_SetItem(argtuple, 0, self_obj);
  for (Py_ssize_t i = 0; i < paramsSize; ++i) {
    PyObject* obj = PyTuple_GET_ITEM(args, i);
    Py_XINCREF(obj);
    PyTuple_SET_ITEM(argtuple, i + 1, obj);
  }
  PyObject* kwdDict = PyDict_Copy(kwds);
  PyDict_SetItemString(kwdDict, "name", self->param_names);
  PyObject* result = PyObject_Call(minuit_function, argtuple, kwdDict);
  Py_XDECREF(argtuple);
  return result;
}

static PyObject *PyFitVCS_enable_energy_spread(PyFitVCSObject *self) {
  self->energy_spread = true;
  return PyLong_FromSsize_t(0);
}

static PyObject *PyFitVCS_disable_energy_spread(PyFitVCSObject *self) {
  self->energy_spread = false;
  return PyLong_FromSsize_t(0);
}

static PyMethodDef PyFitVCS_methods[] = {
  {"minuit", (PyCFunction) PyFitVCS_minuit, METH_VARARGS | METH_KEYWORDS, "Create minuit object"},
  {"enable_energy_spread", (PyCFunction) PyFitVCS_enable_energy_spread, METH_NOARGS, "Enable energy spread"},
  {"disable_energy_spread", (PyCFunction) PyFitVCS_disable_energy_spread, METH_NOARGS, "Disable energy spread"},
  {NULL}
};

static PyTypeObject PyFitVCSType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.FitVCS", /* tp_name */
    sizeof(PyFitVCSObject),  /* tp_basicsize */
    0, /* tp_itemsize */
    (destructor) PyFitVCS_dealloc, /* .tp_dealloc  */
    0, /* tp_print */
    0, /* tp_getattr */
    0, /* tp_setattr */
    0, /* tp_reserved */
    0, /* tp_repr */
    0, /* tp_as_number */
    0, /* tp_as_sequence */
    0, /* tp_as_mapping */
    0, /* tp_hash  */
    PyFitVCS_call, /* tp_call */
    0, /* tp_str */
    0, /* tp_getattro */
    0, /* tp_setattro */
    0, /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "PyFitVCS objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    PyFitVCS_methods, /* tp_methods */
    PyFitVCS_members, /* tp_members */
    0, /* tp_getset  */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) PyFitVCS_init, /* tp_init */
    0, /* tp_alloc */
    PyFitVCS_new, /* tp_new */
};

#endif
