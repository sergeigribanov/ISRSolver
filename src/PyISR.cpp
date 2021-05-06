#include "PyKuraevFadin.hpp"

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
  return PyModule_Create(&PyISR);
}
