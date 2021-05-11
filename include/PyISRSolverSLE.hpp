#ifndef _PY_ISRSOLVER_SLE_HPP_
#define _PY_ISRSOLVER_SLE_HPP_
#include "PyISRSolver.hpp"

static int
PyISRSolverSLE_init(PyISRSolverObject *self, PyObject *args, PyObject *kwds)
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

static PyGetSetDef PyISRSolverSLE_getsetters[] = {
    {"n", (getter) PyISRSolver_n, NULL,
     "Number of points", NULL},
    {NULL}  /* Sentinel */
};

static PyObject *PyISRSolverSLE_bcs_cov_matrix(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[2];
  dims[0] = n;
  dims[1] = n;
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(
      2, dims, NPY_FLOAT64, extractBCSCovMatrix(solver)));
  return PyArray_Transpose(array, NULL);
}

static PyObject *PyISRSolverSLE_bcs_inv_cov_matrix(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[2];
  dims[0] = n;
  dims[1] = n;
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  Eigen::Map<Eigen::MatrixXd> cov(extractBCSCovMatrix(solver), n, n);
  Eigen::MatrixXd icov = cov.inverse();
  const int n2 = n * n;
  double* data = new double[n2];
  std::copy(icov.data(), icov.data() + n2, data);
  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(
      2, dims, NPY_FLOAT64, data));
  return PyArray_Transpose(array, NULL);
}

static PyObject *PyISRSolverSLE_intop_matrix(PyISRSolverObject *self) {
  const std::size_t n = self->solver->getN();
  npy_intp dims[2];
  dims[0] = n;
  dims[1] = n;
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(
      2, dims, NPY_FLOAT64, extractIntOpMatrix(solver)));
  return PyArray_Transpose(array, NULL);
}

static PyMethodDef PyISRSolverSLE_methods[] = {
    {"solve", (PyCFunction) PyISRSolver_solve, METH_NOARGS,
     "Find solution"
    },
    {"save", (PyCFunction) PyISRSolver_save, METH_VARARGS | METH_KEYWORDS,
     "Save results"
    },
    {"bcs", (PyCFunction) PyISRSolver_bcs, METH_NOARGS, "Get numerical solution (Born cross section)"},
    {"ecm", (PyCFunction) PyISRSolver_ecm, METH_NOARGS, "Get center-of-mass energy"},
    {"ecm_err", (PyCFunction) PyISRSolver_ecm_err, METH_NOARGS, "Get center-of-mass energy errors"},
    {"vcs", (PyCFunction) PyISRSolver_vcs, METH_NOARGS, "Get visible cross section"},
    {"vcs_err", (PyCFunction) PyISRSolver_vcs_err, METH_NOARGS, "Get visible cross section errors"},
    {"bcs_cov_matrix", (PyCFunction) PyISRSolverSLE_bcs_cov_matrix, METH_NOARGS,
     "Get covariance matrix of the numerical solution (Born cross section)"},
    {"bcs_inv_cov_matrix", (PyCFunction) PyISRSolverSLE_bcs_inv_cov_matrix, METH_NOARGS,
     "Evaluate inverse covariance matrix of the numerical solution (Born cross section)"},
    {"intop_matrix", (PyCFunction) PyISRSolverSLE_intop_matrix, METH_NOARGS,
     "Integral operator matrix"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyISRSolverSLEType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.PyISRSolverSLE", /* tp_name */
    sizeof(PyISRSolverObject),  /* tp_basicsize */
    0, /* tp_itemsize */
    (destructor) PyISRSolver_dealloc, /* .tp_dealloc  */
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
    PyISRSolverSLE_getsetters, /* tp_getset  */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) PyISRSolverSLE_init, /* tp_init */
    0, /* tp_alloc */
    PyISRSolver_new, /* tp_new */
};

#endif
