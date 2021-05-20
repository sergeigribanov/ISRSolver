#ifndef _PY_FIT_VCS_HPP_
#define _PY_FIT_VCS_HPP_
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <functional>
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>

#include <iostream>

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
  PyObject* params;
  std::function<double(double)> bcsModelLambda;
  std::function<double(double, double)> effLambda;
} PyFitVCSObject;

static void
PyFitVCS_dealloc(PyFitVCSObject *self)
{
  Py_XDECREF(self->energy);
  Py_XDECREF(self->energyErr);
  Py_XDECREF(self->vcs);
  Py_XDECREF(self->vcsErr);
  Py_XDECREF(self->bcsModelFCN);
  Py_XDECREF(self->effFCN);
  Py_XDECREF(self->params);
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
    self->params = NULL;
    self->effLambda =
        [self](double x, double en) {
          PyObject *arglist = Py_BuildValue("(dd)", x, en);
          PyObject *rv = PyObject_CallObject(self->effFCN, arglist);
          double result = PyFloat_AS_DOUBLE(rv);
          Py_CLEAR(rv);
          Py_CLEAR(arglist);
          return result;
        };
    self->bcsModelLambda =
        [self](double en) {
          Py_ssize_t paramsSize = PyTuple_Size(self->params);
          PyObject* argtuple = PyTuple_New(paramsSize + 1);
          PyTuple_SET_ITEM(argtuple, 0, PyFloat_FromDouble(en));
          for (Py_ssize_t i = 0; i < paramsSize; ++i) {
            PyObject* obj = PyTuple_GET_ITEM(self->params, i);
            PyTuple_SET_ITEM(argtuple, i + 1, obj);
          }
          PyObject* rv = PyObject_Call(self->bcsModelFCN, argtuple, NULL);
          double result = PyFloat_AS_DOUBLE(rv);
          Py_XDECREF(rv);
          Py_XDECREF(argtuple);
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
  Py_XINCREF(self->energy);
  Py_XINCREF(self->vcs);
  Py_XINCREF(self->energyErr);
  Py_XINCREF(self->vcsErr);
  Py_XINCREF(self->bcsModelFCN);
  Py_XINCREF(self->effFCN);
  return 0;
}

static PyObject *PyFitVCS_call(PyObject *callable, PyObject *args, PyObject* = NULL) {
  auto self = reinterpret_cast<PyFitVCSObject*>(callable);
  Py_XDECREF(self->params);
  self->params = args;
  Py_XINCREF(self->params);
  // // double tmpResult = self->bcsModelLambda(2.);
  double* energyC = (double*) PyArray_DATA(self->energy);
  // // double* energyErrC = (double*) PyArray_DATA(self->energyErr);
  double* vcsC = (double*) PyArray_DATA(self->vcs);
  // double* vcsErrC = (double*) PyArray_DATA(self->vcsErr);
  // // !!! dimension check
  npy_intp *dims = PyArray_DIMS(self->energy);
  npy_intp dim = dims[0];
  double result = 0;
  for (npy_intp i = 0; i < dim; ++i) {
    double sT = self->threshold * self->threshold;
    double energy = energyC[i];
    std::cout << "energy = " << energy << std::endl;
    std::cout << self->bcsModelLambda(energy) << std::endl;
    // double modelVCS = convolutionKuraevFadin(
   //     energy, self->bcsModelLambda,
   //     0, 1. - sT / energy / energy,
   //     self->effLambda);
  //   double dchi2 = (vcsC[i] - modelVCS) / vcsErrC[i];
  //   dchi2 *= dchi2;
  //   result += dchi2;
  }
  result = 1;
  return PyFloat_FromDouble(result);
}

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
    0, /* tp_methods */
    0, /* tp_members */
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
