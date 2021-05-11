#include "PyKuraevFadin.hpp"
#include "PyISRSolverSLE.hpp"
#include "PyIterISRInterpSolver.hpp"
#include "PyISRSolverTSVD.hpp"
#include "PyISRSolverTikhonov.hpp"

//
//

typedef struct {
    PyObject_HEAD
    PyObject *first; /* first name */
    PyObject *last;  /* last name */
    int number;
} CustomObject;

static void
Custom_dealloc(CustomObject *self)
{
    Py_XDECREF(self->first);
    Py_XDECREF(self->last);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Custom_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    CustomObject *self;
    self = (CustomObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->first = PyUnicode_FromString("");
        if (self->first == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->last = PyUnicode_FromString("");
        if (self->last == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->number = 0;
    }
    return (PyObject *) self;
}

static int
Custom_init(CustomObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = {
    (char*) "first",
    (char*) "last",
    (char*) "number", NULL};
  PyObject *first = NULL, *last = NULL, *tmp;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|UUi", kwlist,
                                   &first, &last,
                                   &self->number))
    return -1;
  if (first) {
    tmp = self->first;
    Py_INCREF(first);
        self->first = first;
        Py_DECREF(tmp);
    }
    if (last) {
        tmp = self->last;
        Py_INCREF(last);
        self->last = last;
        Py_DECREF(tmp);
    }
    return 0;
}

static PyMemberDef Custom_members[] = {
    {"number", T_INT, offsetof(CustomObject, number), 0,
     "custom number"},
    {NULL}  /* Sentinel */
};

static PyObject *
Custom_getfirst(CustomObject *self, void *closure)
{
    Py_INCREF(self->first);
    return self->first;
}

static int
Custom_setfirst(CustomObject *self, PyObject *value, void *closure)
{
    PyObject *tmp;
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the first attribute");
        return -1;
    }
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError,
                        "The first attribute value must be a string");
        return -1;
    }
    tmp = self->first;
    Py_INCREF(value);
    self->first = value;
    Py_DECREF(tmp);
    return 0;
}

static PyObject *
Custom_getlast(CustomObject *self, void *closure)
{
    Py_INCREF(self->last);
    return self->last;
}

static int
Custom_setlast(CustomObject *self, PyObject *value, void *closure)
{
    PyObject *tmp;
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the last attribute");
        return -1;
    }
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError,
                        "The last attribute value must be a string");
        return -1;
    }
    tmp = self->last;
    Py_INCREF(value);
    self->last = value;
    Py_DECREF(tmp);
    return 0;
}

static PyGetSetDef Custom_getsetters[] = {
    {"first", (getter) Custom_getfirst, (setter) Custom_setfirst,
     "first name", NULL},
    {"last", (getter) Custom_getlast, (setter) Custom_setlast,
     "last name", NULL},
    {NULL}  /* Sentinel */
};

static PyObject *
Custom_name(CustomObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("%S %S", self->first, self->last);
}

static PyMethodDef Custom_methods[] = {
    {"name", (PyCFunction) Custom_name, METH_NOARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject CustomType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PyISR.Custom", /* tp_name */
    sizeof(CustomObject),  /* tp_basicsize */
    0, /* tp_itemsize */
    (destructor) Custom_dealloc, /* .tp_dealloc  */
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
    "Custom objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    Custom_methods, /* tp_methods */
    Custom_members, /* tp_members */
    Custom_getsetters, /* tp_getset  */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc) Custom_init, /* tp_init */
    0, /* tp_alloc */
    Custom_new, /* tp_new */
};

//
//

static PyMethodDef PyISR_funcs[] = {
  {"kernelKuraevFadin", (PyCFunction)pyKernelKuraevFadin,
   METH_VARARGS, kernelKuraevFadin_docs},
  {"convolutionKuraevFadin", (PyCFunction)pyConvolutionKuraevFadin,
   METH_VARARGS | METH_KEYWORDS, convolutionKuraevFadin_docs},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef PyISR =
{
  PyModuleDef_HEAD_INIT,
  "PyISR", /* name of module */
  "usage: PyISR.convolutionKuraevFadin(energy, fcn, min_x, max_x)\n", /* module documentation, may be NULL */
  -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
  PyISR_funcs
};

PyMODINIT_FUNC PyInit_PyISR(void)
{
  PyObject *m;
  if (PyType_Ready(&CustomType) < 0)
    return NULL;

  if (PyType_Ready(&PyISRSolverSLEType) < 0)
    return NULL;

  if (PyType_Ready(&PyIterISRInterpSolverType) < 0)
    return NULL;

  if (PyType_Ready(&PyISRSolverTSVDType) < 0)
    return NULL;

  if (PyType_Ready(&PyISRSolverTikhonovType) < 0)
    return NULL;

  m = PyModule_Create(&PyISR);
  if (m == NULL)
    return NULL;

  Py_INCREF(&CustomType);
  if (PyModule_AddObject(m, "Custom", (PyObject *) &CustomType) < 0) {
    Py_DECREF(&CustomType);
    Py_DECREF(m);
    return NULL;
  }

  Py_INCREF(&PyISRSolverSLEType);
  if (PyModule_AddObject(m, "ISRSolverSLE", (PyObject *) &PyISRSolverSLEType) < 0) {
    Py_DECREF(&PyISRSolverSLEType);
    Py_DECREF(m);
    return NULL;
  }

  Py_INCREF(&PyIterISRInterpSolverType);
  if (PyModule_AddObject(m, "IterISRInterpSolver", (PyObject *) &PyIterISRInterpSolverType) < 0) {
    Py_DECREF(&PyIterISRInterpSolverType);
    Py_DECREF(m);
    return NULL;
  }

  Py_INCREF(&PyISRSolverTSVDType);
  if (PyModule_AddObject(m, "ISRSolverTSVD", (PyObject *) &PyISRSolverTSVDType) < 0) {
    Py_DECREF(&PyISRSolverTSVDType);
    Py_DECREF(m);
    return NULL;
  }

  Py_INCREF(&PyISRSolverTikhonovType);
  if (PyModule_AddObject(m, "ISRSolverTikhonov", (PyObject *) &PyISRSolverTikhonovType) < 0) {
    Py_DECREF(&PyISRSolverTikhonovType);
    Py_DECREF(m);
    return NULL;
  }

  import_array();
  return m;
}
