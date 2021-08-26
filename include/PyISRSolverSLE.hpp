#ifndef _PY_ISRSOLVER_SLE_HPP_
#define _PY_ISRSOLVER_SLE_HPP_
#include "Chi2Test.hpp"
#include "PyISRSolver.hpp"

static PyObject *PyISRSolverSLE_set_interp_settings(PyISRSolverObject *self, PyObject *args) {
  // !!! TO-DO: return none
  PyObject* obj;
  if (!PyArg_ParseTuple(args, "O", &obj)) {
    return 0;
  }
  PyObject *iterator = PyObject_GetIter(obj);
  PyObject* item;
  PyObject* subitem;
  PyObject* index;
  bool el0;
  int el1;
  int el2;
  std::vector<std::tuple<bool, int, int>> result;
  result.reserve(PyObject_Length(obj));
  while ((item = PyIter_Next(iterator))) {
    index = PyLong_FromSsize_t(0);
    subitem = PyObject_GetItem(item, index);
    if (PyObject_IsTrue(subitem)) {
      el0 = true;
    } else {
      el0 = false;
    }
    Py_DECREF(index);
    Py_DECREF(subitem);

    index = PyLong_FromSsize_t(1);
    subitem = PyObject_GetItem(item, index);
    el1 = PyLong_AsLong(subitem);
    Py_DECREF(index);
    Py_DECREF(subitem);

    index = PyLong_FromSsize_t(2);
    subitem = PyObject_GetItem(item, index);
    el2 = PyLong_AsLong(subitem);
    Py_DECREF(index);
    Py_DECREF(subitem);
    result.push_back(std::make_tuple(el0, el1, el2));
    Py_DECREF(item);
  }
  Py_DECREF(iterator);
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  solver->setRangeInterpSettings(result);
  return PyLong_FromSsize_t(0);
}

static PyGetSetDef PyISRSolverSLE_getsetters[] = {
    {"n", (getter) PyISRSolver_n, NULL,
     "Number of points", NULL},
    {"energy_spread_enabled", (getter) PyISRSolver_get_energy_sread_enabled,
     (setter) PyISRSolverTikhonov_set_energy_sread_enabled, "Energy spread flag", NULL},
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

static PyObject *PyISRSolverSLE_interp_eval(PyISRSolverObject *self, PyObject *args, PyObject *kwds) {
  double energy;
  static const char *kwlist[] = {"cm_energy", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d",
                                   const_cast<char**>(kwlist),
                                   &energy)) {
    return 0;
  }
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  const double value = solver->interpEval(solver->bcs(), energy);
  return PyFloat_FromDouble(value);
}


static PyObject *PyISRSolverSLE_condnum_eval(PyISRSolverObject *self) {
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  const double value = solver->evalConditionNumber();
  return PyFloat_FromDouble(value);
}

static PyObject *PyISRSolverSLE_eval_equation_matrix(PyISRSolverObject *self) {
  ISRSolverSLE* solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  solver->evalEqMatrix();
  return PyLong_FromSsize_t(0);
}

static PyObject *PyISRSolverSLE_chi2_test_model(PyISRSolverObject *self, PyObject *args, PyObject *kwds) {
  int n;
  double ampl;
  char *modelPath = NULL;
  char *modelVCsName = NULL;
  char *modelBCsName = NULL;
  char *outputPath = NULL;
  static const char *kwlist[] =
    {"n", "initial_ampl", "model_path", "model_visible_cs_name",
    "model_born_cs_name", "output_path", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "idssss",
                                   const_cast<char**>(kwlist),
                                   &n, &ampl,
                                   &modelPath,
                                   &modelVCsName,
                                   &modelBCsName,
                                   &outputPath)) {
    return 0;
  }
  auto solver = reinterpret_cast<ISRSolverSLE*>(self->solver);
  chi2TestModel(solver,
                {.n = n,
                 .initialChi2Ampl = ampl,
                 .modelPath = modelPath,
                 .modelVCSName = modelVCsName,
                 .modelBCSName = modelBCsName,
                 .outputPath = outputPath});
  return PyLong_FromSsize_t(0);
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
    {"condnum_eval", (PyCFunction) PyISRSolverSLE_condnum_eval, METH_NOARGS,
     "Eval condition number of the integral operator matrix F or GF in the case, when c.m. energy spread is enabled."},
    {"interp_eval", (PyCFunction) PyISRSolverSLE_interp_eval, METH_VARARGS | METH_KEYWORDS,
     "Eval interpolation value at certain c.m. energy point"},
    {"eval_equation_matrix", (PyCFunction) PyISRSolverSLE_eval_equation_matrix, METH_NOARGS,
    "Eval equation matrix"},
    {"reset_ecm_err", (PyCFunction) PyISRSolver_reset_ecm_err,  METH_VARARGS | METH_KEYWORDS,
     "Reset c.m. energy spread"},
    {"set_interp_settings", (PyCFunction) PyISRSolverSLE_set_interp_settings, METH_VARARGS, "Set interpolation settings"},
    {"chi2_test_model", (PyCFunction) PyISRSolverSLE_chi2_test_model, METH_VARARGS | METH_KEYWORDS,
    "Chi2 test using model Born and Visible cross sections"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyISRSolverSLEType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.ISRSolverSLE", /* tp_name */
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
    "ISRSolverSLE objects", /* tp_doc */
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
    (initproc) PyISRSolver_init<ISRSolverSLE>, /* tp_init */
    0, /* tp_alloc */
    PyISRSolver_new, /* tp_new */
};

#endif
