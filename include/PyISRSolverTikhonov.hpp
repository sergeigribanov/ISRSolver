#ifndef _PY_ISRSOLVER_TIKHONOV_HPP_
#define _PY_ISRSOLVER_TIKHONOV_HPP_
#include <vector>
#include <nlopt.hpp>
#include "PyISRSolverSLE.hpp"
#include "ISRSolverTikhonov.hpp"
#include "Utils.hpp"

#include <iostream>

static PyObject *
PyISRSolverTikhonov_get_lambda(PyISRSolverObject *self, void *closure)
{
  auto solver = reinterpret_cast<ISRSolverTikhonov*>(self->solver);
  return PyFloat_FromDouble(solver->getLambda());
}

static int
PyISRSolverTikhonov_set_lambda(PyISRSolverObject *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete regularization parameter");
    return -1;
  }
  if (!PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError,
                    "Regularization parameter must be float");
    return -1;
  }
  auto solver = reinterpret_cast<ISRSolverTikhonov*>(self->solver);
  solver->setLambda(PyFloat_AS_DOUBLE(value));
  return 0;
}

static PyGetSetDef PyISRSolverTikhonov_getsetters[] = {
  {"n", (getter) PyISRSolver_n, NULL, "Number of points", NULL},
  {"energy_spread_enabled", (getter) PyISRSolver_get_energy_sread_enabled,
   (setter) PyISRSolverTikhonov_set_energy_sread_enabled, "Energy spread flag", NULL},
  {"reg_param", (getter) PyISRSolverTikhonov_get_lambda,
   (setter) PyISRSolverTikhonov_set_lambda, "Regularization parameter", NULL},
  {NULL}  /* Sentinel */
};

static PyObject *PyISRSolverTikhonov_make_LCurve(PyISRSolverObject *self, PyObject *args, PyObject *kwds) {
  PyArrayObject *lambdas = NULL;
  static const char *kwlist[] = {"lambdas", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!",
                                   const_cast<char**>(kwlist),
                                   &PyArray_Type, &lambdas)) {
    return 0;
  }
  const int ndim = PyArray_NDIM(lambdas);
  if (ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Wrong dimension of the array with regularization parameter values. The array must be 1D array.");
    return 0;
  }
  npy_intp *dims = PyArray_DIMS(lambdas);
  npy_intp dim = dims[0];
  double *lambdasC = (double*) PyArray_DATA(lambdas);
  auto solver = reinterpret_cast<ISRSolverTikhonov*>(self->solver);
  double* x = new double[dim];
  double* y = new double[dim];
  double* curv = new double[dim];
  for (npy_intp i = 0; i < dim; ++i) {
    solver->setLambda(lambdasC[i]);
    solver->solve();
    x[i] = solver->evalEqNorm2();
    y[i] = solver->evalSmoothnessConstraintNorm2();
    curv[i] = -solver->evalLCurveCurvature();
  }
  PyObject* aX = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, x);
  PyObject* aY = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, y);
  PyObject* aC = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, curv);
  PyObject* result = PyDict_New();
  PyDict_SetItemString(result, "x", aX);
  PyDict_SetItemString(result, "y", aY);
  PyDict_SetItemString(result, "c", aC);
  return result;
}

static PyObject *PyISRSolverTikhonov_solve_LCurve(PyISRSolverObject *self, PyObject *args, PyObject *kwds) {
  double lambda;
  static const char *kwlist[] = {"initial_lambda", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d",
                                   const_cast<char**>(kwlist),
                                   &PyArray_Type, &lambda)) {
    return 0;
  }
  /**
   * Creating NLOPT optimizer
   */
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> z(1, lambda);
  std::vector<double> lowerBounds(1, 0);
  opt.set_lower_bounds(lowerBounds);
  opt.set_min_objective(lambdaObjective, self->solver);
  opt.set_xtol_rel(1.e-6);
  double minf;
  /**
   * Run the L-curve curvature maximization
   */
  opt.optimize(z, minf);
  PyObject* pyLambda = PyFloat_FromDouble(z[0]);
  PyObject* pyMaxCurv = PyFloat_FromDouble(-minf);
  PyObject* result = PyDict_New();
  PyDict_SetItemString(result, "lambda", pyLambda);
  PyDict_SetItemString(result, "max_curv", pyMaxCurv);
  return result;
}

static PyMethodDef PyISRSolverTikhonov_methods[] = {
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
    {"make_LCurve", (PyCFunction) PyISRSolverTikhonov_make_LCurve, METH_VARARGS | METH_KEYWORDS,
     "Obtaining L-curve and it's curvature"},
    {"solve_LCurve", (PyCFunction) PyISRSolverTikhonov_solve_LCurve, METH_VARARGS | METH_KEYWORDS,
    "Find solution using the L-curve criterion"},
    {"chi2_test_model", (PyCFunction) PyISRSolverSLE_chi2_test_model, METH_VARARGS | METH_KEYWORDS,
     "Chi2 test of the numerical solution with respect to the modle Born cross section. Visible cross section is generated multiple times using model visible cross section."},
    {"ratio_test_model", (PyCFunction) PyISRSolverSLE_ratio_test_model, METH_VARARGS | METH_KEYWORDS,
     "Ratio test (test of interpolation quality). Ratio of the numerical solution to the model Born cross section. Visible cross section is generated multiple times using model visible cross section."},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyISRSolverTikhonovType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.ISRSolverTikhonov", /* tp_name */
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
    "ISRSolverTikhonov objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    PyISRSolverTikhonov_methods, /* tp_methods */
    0, //PyISRSolverSLE_members, /* tp_members */
    PyISRSolverTikhonov_getsetters, /* tp_getset  */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) PyISRSolver_init<ISRSolverTikhonov>, /* tp_init */
    0, /* tp_alloc */
    PyISRSolver_new, /* tp_new */
};

#endif
