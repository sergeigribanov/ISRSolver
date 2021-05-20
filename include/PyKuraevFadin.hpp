#ifndef _PY_KURAEV_FADIN_HPP_
#define _PY_KURAEV_FADIN_HPP_
#define PY_SSIZE_T_CLEAN
#include <functional>
#include <Python.h>
#include <structmember.h>
#include "KuraevFadin.hpp"
static PyObject* pyKernelKuraevFadin(PyObject* self,
                                   PyObject* args) {
  double xC;
  double sC;
  if (!PyArg_ParseTuple(args, "dd", &xC, &sC)) {
    return 0;
  }
  return PyFloat_FromDouble(kernelKuraevFadin(xC, sC));
}

static PyObject* pyConvolutionKuraevFadin(PyObject* self,
                                          PyObject* args,
                                          PyObject* kwds) {
  PyObject *cb = nullptr;
  PyObject *efficiency = nullptr;
  double energyC;
  double minXC;
  double maxXC;
  static const char *kwlist[] = {"energy", "fcn", "min_x", "max_x", "efficiency", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "dOdd|O",
                                   const_cast<char**>(kwlist),
                                   &energyC, &cb, &minXC, &maxXC, &efficiency)) {
    return 0;
  }
  if (!PyCallable_Check(cb)) {
    PyErr_SetString(PyExc_TypeError, "convolutionKuraevFadin: a callable is required");
    return 0;
  }
  std::function<double(double)> fcnC =
      [cb](double en) {
        PyObject *arglist = Py_BuildValue("(d)", en);;
        PyObject *rv = PyObject_CallObject(cb, arglist);
        double result = PyFloat_AS_DOUBLE(rv);
        Py_CLEAR(rv);
        Py_CLEAR(arglist);
        return result;
      };
  if (!efficiency) {
    double result = convolutionKuraevFadin(energyC, fcnC, minXC, maxXC);
    return PyFloat_FromDouble(result);
  }
  std::function<double(double, double)> efficiencyC =
      [efficiency](double x, double en) {
        PyObject *arglist = Py_BuildValue("(dd)", x, en);;
        PyObject *rv = PyObject_CallObject(efficiency, arglist);
        double result = PyFloat_AS_DOUBLE(rv);
        Py_CLEAR(rv);
        Py_CLEAR(arglist);
        return result;
      };
  double result = convolutionKuraevFadin(energyC, fcnC, minXC, maxXC, efficiencyC);
  return PyFloat_FromDouble(result);
}

static char kernelKuraevFadin_docs[] =
    "kernelKuraevFadin(x, s): Any message you want to put here!!\n";

static char convolutionKuraevFadin_docs[] =
    "convolutionKuraevFadin(energy, fcn, min_x, max_x, efficiency = lambda x, energy: 1.): Any message you want to put here!!\n";

#endif
