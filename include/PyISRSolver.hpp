#ifndef _PY_ISRSOLVER_HPP_
#define _PY_ISRSOLVER_HPP_
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include "ISRSolverSLE.hpp"

typedef struct {
  PyObject_HEAD
  PyObject *eff;
  BaseISRSolver *solver;
} PyISRSolverObject;

static void
PyISRSolver_dealloc(PyISRSolverObject *self)
{
  if (self->solver) {
    delete self->solver;
  }
  Py_XDECREF(self->eff);
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
PyISRSolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyISRSolverObject *self;
  self = reinterpret_cast<PyISRSolverObject*>(type->tp_alloc(type, 0));
  if (self != NULL) {
    self->solver = nullptr;
  }
  return reinterpret_cast<PyObject*>(self);
}

static PyObject *
PyISRSolver_n(PyISRSolverObject *self, void *closure)
{
  return PyLong_FromSsize_t(self->solver->getN());
}

static PyObject *PyISRSolver_ecm(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[1];
  dims[0] = n;
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, extractECMPointer(self->solver));
}

static PyObject *PyISRSolver_ecm_err(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[1];
  dims[0] = n;
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, extractECMErrPointer(self->solver));
}

static PyObject *PyISRSolver_vcs(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[1];
  dims[0] = n;
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, extractVCSPointer(self->solver));
}

static PyObject *PyISRSolver_vcs_err(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[1];
  dims[0] = n;
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, extractVCSErrPointer(self->solver));
}

static PyObject *PyISRSolver_bcs(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[1];
  dims[0] = n;
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, extractBCSPointer(self->solver));
}

static PyObject *PyISRSolver_solve(PyISRSolverObject *self) {
  // !!! TO-DO: return none
  self->solver->solve();
  return PyLong_FromSsize_t(0);
}

static PyObject *PyISRSolver_save(PyISRSolverObject *self, PyObject *args, PyObject *kwds) {
  char* outputPath = NULL;
  char* visibleCSGraphName = NULL;
  char* bornCSGraphName = NULL;
  static const char *kwlist[] = {"output_path", "vcs_name", "bcs_name", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss",
                                   const_cast<char**>(kwlist),
                                   &outputPath,
                                   &visibleCSGraphName,
                                   &bornCSGraphName)) {
    return 0;
  }
  std::string outputPathS = outputPath;
  std::string visibleCSGraphNameS;
  std::string bornCSGraphNameS;
  //!!!  delete [] outputPath;
  if (!visibleCSGraphName) {
    visibleCSGraphNameS = "vcs";
  } else {
    visibleCSGraphNameS = visibleCSGraphName;
    //!!! delete [] visibleCSGraphName;
  }
  if (!bornCSGraphName) {
    bornCSGraphNameS = "bcs";
  } else {
    bornCSGraphNameS = bornCSGraphName;
    //!!! delete [] bornCSGraphName;
  }
  self->solver->save(outputPathS,
                     {.visibleCSGraphName = visibleCSGraphNameS,
                      .bornCSGraphName = bornCSGraphNameS});
  return PyLong_FromSsize_t(0);
}


#endif
