#ifndef _PY_ISRSOLVER_SLE_HPP_
#define _PY_ISRSOLVER_SLE_HPP_
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include "ISRSolverSLE.hpp"

typedef struct {
    PyObject_HEAD
    PyObject *eff;
    ISRSolverSLE *solver;
} PyISRSolverSLEObject;

static void
PyISRSolverSLE_dealloc(PyISRSolverSLEObject *self)
{
  if (self->solver) {
    delete self->solver;
  }
  Py_XDECREF(self->eff);
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
PyISRSolverSLE_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyISRSolverSLEObject *self;
    self = (PyISRSolverSLEObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
      self->solver = nullptr;
    }
    return (PyObject *) self;
}

static int
PyISRSolverSLE_init(PyISRSolverSLEObject *self, PyObject *args, PyObject *kwds)
{
  // !!! TO-DO: check dimension
  unsigned long nC;
  PyArrayObject *energy = NULL;
  PyArrayObject *visibleCS = NULL;
  PyArrayObject *energyErr = NULL;
  PyArrayObject *visibleCSErr = NULL;
  PyObject* enableEnergySpread = NULL;
  double thresholdC;
  self->eff = NULL;
  static const char *kwlist[] =
      {"n", "energy", "vis_cs", "energy_err",
       "vis_cs_err", "threshold", "efficiency",
       "enableEnergySpread", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "kO!O!O!O!d|OO",
                                   const_cast<char**>(kwlist),
                                   &nC,
                                   &PyArray_Type, &energy,
                                   &PyArray_Type, &visibleCS,
                                   &PyArray_Type, &energyErr,
                                   &PyArray_Type, &visibleCSErr,
                                   &thresholdC, &self->eff,
                                   &enableEnergySpread)) {
    return -1;
  }
  Py_XINCREF(self->eff);
  if (!PyCallable_Check(self->eff)) {
    PyErr_SetString(PyExc_TypeError, "ISRSolverSLE: a callable efficiency object is required");
    return -1;
  }
  std::function<double(double, double)> effC =
      [self](double x, double en) {
        PyObject *arglist = Py_BuildValue("(dd)", x, en);;
        PyObject *rv = PyObject_CallObject(self->eff, arglist);
        double result = PyFloat_AS_DOUBLE(rv);
        Py_CLEAR(rv);
        Py_CLEAR(arglist);
        return result;
      };
  double *energyCArrayC = (double*) PyArray_DATA(energy);
  double *visibleCSArrayC = (double*) PyArray_DATA(visibleCS);
  double *energyErrArrayC = (double*) PyArray_DATA(energyErr);
  double *visibleCSErrArrayC = (double*) PyArray_DATA(visibleCSErr);
  self->solver = new ISRSolverSLE(nC,
                                  energyCArrayC, visibleCSArrayC,
                                  energyErrArrayC, visibleCSErrArrayC,
                                  thresholdC, effC);
  if (!enableEnergySpread) {
    self->solver->disableEnergySpread();
  } else if (PyObject_IsTrue(enableEnergySpread)) {
    self->solver->enableEnergySpread();
  } else {
    self->solver->disableEnergySpread();
  }
  return 0;
}

// static PyMemberDef PyISRSolverSLE_members[] = {
//     {"number", T_INT, offsetof(PyISRSolverSLEObject, number), 0,
//      "custom number"},
//     {NULL}  /* Sentinel */
// };

// static PyObject *
// PyISRSolverSLE_getfirst(PyISRSolverSLEObject *self, void *closure)
// {
//     Py_INCREF(self->first);
//     return self->first;
// }

// static int
// PyISRSolverSLE_setfirst(PyISRSolverSLEObject *self, PyObject *value, void *closure)
// {
//     PyObject *tmp;
//     if (value == NULL) {
//         PyErr_SetString(PyExc_TypeError, "Cannot delete the first attribute");
//         return -1;
//     }
//     if (!PyUnicode_Check(value)) {
//         PyErr_SetString(PyExc_TypeError,
//                         "The first attribute value must be a string");
//         return -1;
//     }
//     tmp = self->first;
//     Py_INCREF(value);
//     self->first = value;
//     Py_DECREF(tmp);
//     return 0;
// }

// static PyObject *
// PyISRSolverSLE_getlast(PyISRSolverSLEObject *self, void *closure)
// {
//     Py_INCREF(self->last);
//     return self->last;
// }

// static int
// PyISRSolverSLE_setlast(PyISRSolverSLEObject *self, PyObject *value, void *closure)
// {
//     PyObject *tmp;
//     if (value == NULL) {
//         PyErr_SetString(PyExc_TypeError, "Cannot delete the last attribute");
//         return -1;
//     }
//     if (!PyUnicode_Check(value)) {
//         PyErr_SetString(PyExc_TypeError,
//                         "The last attribute value must be a string");
//         return -1;
//     }
//     tmp = self->last;
//     Py_INCREF(value);
//     self->last = value;
//     Py_DECREF(tmp);
//     return 0;
// }

// static PyGetSetDef PyISRSolverSLE_getsetters[] = {
//     {"first", (getter) PyISRSolverSLE_getfirst, (setter) PyISRSolverSLE_setfirst,
//      "first name", NULL},
//     {"last", (getter) PyISRSolverSLE_getlast, (setter) PyISRSolverSLE_setlast,
//      "last name", NULL},
//     {NULL}  /* Sentinel */
// };

// static PyObject *
// PyISRSolverSLE_name(PyISRSolverSLEObject *self, PyObject *Py_UNUSED(ignored))
// {
//     return PyUnicode_FromFormat("%S %S", self->first, self->last);
// }

// static PyMethodDef PyISRSolverSLE_methods[] = {
//     {"name", (PyCFunction) PyISRSolverSLE_name, METH_NOARGS,
//      "Return the name, combining the first and last name"
//     },
//     {NULL}  /* Sentinel */
// };


static PyObject *PyISRSolverSLE_solve(PyISRSolverSLEObject *self) {
  // !!! TO-DO: return none
  self->solver->solve();
  return PyLong_FromSsize_t(0);
}

static PyObject *PyISRSolverSLE_save(PyISRSolverSLEObject *self, PyObject *args, PyObject *kwds) {
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

static PyMethodDef PyISRSolverSLE_methods[] = {
    {"solve", (PyCFunction) PyISRSolverSLE_solve, METH_NOARGS,
     "Find solution"
    },
    {"save", (PyCFunction) PyISRSolverSLE_save, METH_VARARGS | METH_KEYWORDS,
     "Save results"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject PyISRSolverSLEType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.PyISRSolverSLE", /* tp_name */
    sizeof(PyISRSolverSLEObject),  /* tp_basicsize */
    0, /* tp_itemsize */
    (destructor) PyISRSolverSLE_dealloc, /* .tp_dealloc  */
    0, /* tp_print */
    0, /* tp_getattr */
    0, /* tp_setattr */
    0, /* tp_reserved */
    0, /* tp_repr */
    0, /* tp_as_number */
    0, /* tp_as_sequence */
    0, /* tp_as_mapping */
    0, /* tp_hash  */
    0, /* tp_call */
    0, /* tp_str */
    0, /* tp_getattro */
    0, /* tp_setattro */
    0, /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "PyISRSolverSLE objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    PyISRSolverSLE_methods, /* tp_methods */
    0, //PyISRSolverSLE_members, /* tp_members */
    0, //PyISRSolverSLE_getsetters, /* tp_getset  */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) PyISRSolverSLE_init, /* tp_init */
    0, /* tp_alloc */
    PyISRSolverSLE_new, /* tp_new */
};

#endif
