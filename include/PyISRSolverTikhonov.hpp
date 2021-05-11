#ifndef _PY_ISRSOLVER_TIKHONOV_HPP_
#define _PY_ISRSOLVER_TIKHONOV_HPP_
#include "PyISRSolverSLE.hpp"
#include "ISRSolverTikhonov.hpp"

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
    {"set_interp_settings", (PyCFunction) PyISRSolverSLE_set_interp_settings, METH_VARARGS, "Set interpolation settings"},
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
