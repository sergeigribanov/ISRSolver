#ifndef _PY_ITER_ISR_SOLVER_USE_VCS_FIT_HPP_
#define _PY_ITER_ISR_SOLVER_USE_VCS_FIT_HPP_
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <Eigen/Dense>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>
#include <functional>
#include "Integration.hpp"
#include "KuraevFadin.hpp"

#include <iostream>

typedef struct {
  PyObject_HEAD
  PyArrayObject* energy;
  PyArrayObject* vcs;
  PyArrayObject* vcsErr;
  PyObject* vcsFitFCN;
  PyObject* effFCN;
  std::function<double(double, double)> effLambda;
  bool energy_spread;
  double sigmaEn;
  double threshold;
  unsigned int npoints;
  unsigned int niter;
  Eigen::VectorXd rad_corr;
  Eigen::VectorXd bcs;
  Eigen::VectorXd bcsErr;
} PyIterISRSolverUseVCSFitObject;

static void
PyIterISRSolverUseVCSFit_dealloc(PyIterISRSolverUseVCSFitObject *self)
{
  Py_XDECREF(self->energy);
  Py_XDECREF(self->vcs);
  Py_XDECREF(self->vcsErr);
  Py_XDECREF(self->vcsFitFCN);
  Py_XDECREF(self->effFCN);
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
PyIterISRSolverUseVCSFit_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyIterISRSolverUseVCSFitObject *self;
  self = reinterpret_cast<PyIterISRSolverUseVCSFitObject*>(type->tp_alloc(type, 0));
  if (self != NULL) {
    self->energy = NULL;
    self->vcs = NULL;
    self->vcsErr = NULL;
    self->vcsFitFCN = NULL;
    self->effFCN = NULL;
    self->npoints = 100;
    self->niter = 10;
    self->energy_spread = false;
    self->sigmaEn = 0.;
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
PyIterISRSolverUseVCSFit_init(PyIterISRSolverUseVCSFitObject *self, PyObject *args, PyObject *kwds) {
  static const char *kwlist[] =
      {"threshold", "energy", "vcs", "vcs_err", "vcs_fcn", "efficiency", "energy_spread", NULL};
  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "dO!O!O!OO|d",
          const_cast<char**>(kwlist),
          &(self->threshold),
          &PyArray_Type, &(self->energy),
          &PyArray_Type, &(self->vcs),
          &PyArray_Type, &(self->vcsErr),
          &(self->vcsFitFCN), &(self->effFCN),
          &(self->sigmaEn))) {
    return -1;
  }
  // TO-DO: check callables !!!
  Py_XINCREF(self->energy);
  Py_XINCREF(self->vcs);
  Py_XINCREF(self->vcsErr);
  Py_XINCREF(self->vcsFitFCN);
  Py_XINCREF(self->effFCN);
  return 0;
}

static PyObject *PyIterISRSolverUseFit_enable_energy_spread(
    PyIterISRSolverUseVCSFitObject *self) {
  self->energy_spread = true;
  return PyLong_FromSsize_t(0);
}

static PyObject *PyIterISRSolverUseFit_disable_energy_spread(
    PyIterISRSolverUseVCSFitObject *self) {
  self->energy_spread = false;
  return PyLong_FromSsize_t(0);
}

static PyObject *PyIterISRSolverUseFit_solve(
    PyIterISRSolverUseVCSFitObject *self, PyObject *args, PyObject *kwds) {
  static const char *kwlist[] =
      {"verbose",  NULL};
  bool verbose = false;
  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|p",
          const_cast<char**>(kwlist),
          &verbose)) {
    return 0;
  }
  Eigen::VectorXd ecm = Eigen::VectorXd::Zero(self->npoints);
  PyObject *maxelem = PyArray_Max(self->energy, 0, NULL);
  const double maxen = PyFloat_AS_DOUBLE(maxelem);
  Py_XDECREF(maxelem);
  const double eh = (maxen - self->threshold) / (self->npoints - 1.);
  for (unsigned int i = 0; i < self->npoints; ++i) {
    ecm(i) = self->threshold + eh * i;
  }
  Eigen::VectorXd radCorr = Eigen::VectorXd::Zero(self->npoints);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, self->npoints);
  gsl_spline_init (spline, ecm.data(), radCorr.data(), self->npoints);
  std::function<double(double)> radFCN =
      [&spline, &acc](double en) {
        double result = gsl_spline_eval(spline, en, acc);
        return result;
      };
  /**
   * Born cross section function
   */
  std::function<double(double)> born_fcn =
      [maxen, &radFCN, &self](double en) {
        if (en <= self->threshold) {
          return 0.;
        }
        double result = 0;
        if (en >= maxen) {
          PyObject* argtuple = PyTuple_New(1);
          PyTuple_SET_ITEM(argtuple, 0, PyFloat_FromDouble(maxen));
          PyObject* rv = PyObject_Call(self->vcsFitFCN, argtuple, NULL);
          const double tmp_val = PyFloat_AS_DOUBLE(rv);
          Py_XDECREF(rv);
          result = tmp_val / (1. + radFCN(maxen));
          return result;
        }
        PyObject* argtuple = PyTuple_New(1);
        PyTuple_SET_ITEM(argtuple, 0, PyFloat_FromDouble(en));
        PyObject* rv = PyObject_Call(self->vcsFitFCN, argtuple, NULL);
        const double tmp_val = PyFloat_AS_DOUBLE(rv);
        Py_XDECREF(rv);
        result = tmp_val / (1. + radFCN(en));
        return result;
      };
  std::function<double(double)> vcs_fcn_no_spread =
      [&self, &born_fcn](double en) {
        if (en <= self->threshold) {
          return 0.;
        }
        const double s_threshold = self->threshold * self->threshold;
        const double s = en * en;
        double result = 0;
        if (self->effFCN) {
          result = convolutionKuraevFadin(en, born_fcn, 0, 1. - s_threshold / s,
                                          self->effLambda);
        } else {
          result = convolutionKuraevFadin(en, born_fcn, 0, 1. - s_threshold / s);
        }
        return result;
      };
  std::function<double(double)> vcs_fcn =
      [&self, &vcs_fcn_no_spread](double en) {
        double result = 0;
        if (self->energy_spread) {
          result = gaussian_conv(en, self->sigmaEn * self->sigmaEn,
                                 vcs_fcn_no_spread);
        } else {
          result = vcs_fcn_no_spread(en);
        }
        return result;
      };
  Eigen::VectorXd tmpRad = Eigen::VectorXd::Zero(self->npoints);
  for (unsigned iter = 0; iter < self->niter; ++iter) {
    if (verbose) {
      std::cout << "ITER: " << iter << " / " << self->niter << std::endl;
    }
    for (unsigned int i = 0; i < self->npoints; ++i) {
      tmpRad(i) = vcs_fcn(ecm(i)) / born_fcn(ecm(i)) - 1;
      if (std::isnan(tmpRad(i))) {
        tmpRad(i) = 0;
      }
    }
    std::copy(tmpRad.data(), tmpRad.data() + self->npoints, radCorr.data());
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_linear, self->npoints);
    gsl_spline_init(spline, ecm.data(), radCorr.data(), self->npoints);
  }
  npy_intp *dims = PyArray_DIMS(self->energy);
  npy_intp dim = dims[0];
  self->rad_corr = Eigen::VectorXd(dim);
  self->bcs = Eigen::VectorXd(dim);
  self->bcsErr = Eigen::VectorXd(dim);
  double* energyC = (double*) PyArray_DATA(self->energy);
  double* vcsC = (double*) PyArray_DATA(self->vcs);
  double* vcsErrC = (double*) PyArray_DATA(self->vcsErr);
  for (npy_intp i = 0; i < dim; ++i) {
    double energy = energyC[i];
    (self->rad_corr)(i) = radFCN(energy);
    (self->bcs)(i) = vcsC[i] / (1. + (self->rad_corr)(i));
    (self->bcsErr)(i) = vcsErrC[i] / (1. + (self->rad_corr)(i));
  }
  return PyLong_FromSsize_t(0);
}

static PyObject *PyIterISRSolverUseFit_rad_corr(PyIterISRSolverUseVCSFitObject *self) {
  npy_intp dims[1];
  dims[0] = self->rad_corr.size();
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, self->rad_corr.data());
}

static PyObject *PyIterISRSolverUseFit_ecm(PyIterISRSolverUseVCSFitObject *self) {
  PyObject* obj = reinterpret_cast<PyObject*>(self->energy);
  Py_XINCREF(obj);
  return obj;
}

static PyObject *PyIterISRSolverUseFit_vcs(PyIterISRSolverUseVCSFitObject *self) {
  PyObject* obj = reinterpret_cast<PyObject*>(self->vcs);
  Py_XINCREF(obj);
  return obj;
}

static PyObject *PyIterISRSolverUseFit_vcs_err(PyIterISRSolverUseVCSFitObject *self) {
  PyObject* obj = reinterpret_cast<PyObject*>(self->vcsErr);
  Py_XINCREF(obj);
  return obj;
}

static PyObject *PyIterISRSolverUseFit_bcs(PyIterISRSolverUseVCSFitObject *self) {
  npy_intp dims[1];
  dims[0] = self->bcs.size();
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, self->bcs.data());
}

static PyObject *PyIterISRSolverUseFit_bcs_err(PyIterISRSolverUseVCSFitObject *self) {
  npy_intp dims[1];
  dims[0] = self->bcsErr.size();
  return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, self->bcsErr.data());
}

static PyMethodDef PyIterISRSolverUseFit_methods[] = {
  {"enable_energy_spread", (PyCFunction) PyIterISRSolverUseFit_enable_energy_spread,
   METH_NOARGS, "Enable energy spread"},
  {"disable_energy_spread", (PyCFunction) PyIterISRSolverUseFit_disable_energy_spread,
   METH_NOARGS, "Disable energy spread"},
  {"solve", (PyCFunction) PyIterISRSolverUseFit_solve,
   METH_VARARGS | METH_KEYWORDS, "Find solution"},
  {"ecm", (PyCFunction) PyIterISRSolverUseFit_ecm, METH_NOARGS, "Get center-of-mass energy"},
  {"vcs", (PyCFunction) PyIterISRSolverUseFit_vcs, METH_NOARGS, "Get visible cross section"},
  {"vcs_err", (PyCFunction) PyIterISRSolverUseFit_vcs_err, METH_NOARGS, "Get visible cross section error"},
  {"rad_corr", (PyCFunction) PyIterISRSolverUseFit_rad_corr, METH_NOARGS, "Get radiative correction"},
  {"bcs", (PyCFunction) PyIterISRSolverUseFit_bcs, METH_NOARGS, "Get Born cross section"},
  {"bcs_err", (PyCFunction) PyIterISRSolverUseFit_bcs_err, METH_NOARGS, "Get Born cross section error"},
  {NULL}
};

static PyTypeObject PyIterISRSolverUseVCSFitType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "PyISR.IterISRSolverUseVCSFit", /* tp_name */
  sizeof(PyIterISRSolverUseVCSFitObject), /* tp_basicsize */
  0, /* tp_itemsize */
  (destructor) PyIterISRSolverUseVCSFit_dealloc, /* tp_dealloc */
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
  "PyIterISRSolverUseVCSFit objects", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  PyIterISRSolverUseFit_methods, /* tp_methods */
  0, /* tp_members */
  0, /* tp_getset  */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) PyIterISRSolverUseVCSFit_init, /* tp_init */
  0, /* tp_alloc */
  PyIterISRSolverUseVCSFit_new, /* tp_new */
};

#endif
