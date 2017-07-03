#include <itkMesh.h>
#include <itkTransformMeshFilter.h>
#include <itkCompositeTransform.h>
#include <itkLBFGSBOptimizer.h>
#include <itkImageRegistrationMethodv4.h>

#include "itkGMMPointSetToPointSetRegistrationMethod.h"
#include "itkPointSetPropertiesCalculator.h"
#include "itkNormalizePointSet.h"
#include "itkInitializeTransform.h"
#include "itkInitializeMetric.h"
#include "itkPointSetToPointSetMetrics.h"

#include "args.hxx"
#include "argsCustomParsers.h"

#include "agtkIO.h"
#include "agtkCommandIterationUpdate.h"

using namespace agtk;

typedef itk::Mesh<float, 3U> MeshType;
typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
typedef itk::Transform <double, MeshType::PointDimension, MeshType::PointDimension> TransformType;

int main(int argc, char** argv) {
  args::ArgumentParser parser("GMM PointSet Registration", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Group allRequired(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argFixedFileName(allRequired, "fixed", "The fixed mesh (point-set) filename", {'f', "fixed"});
  args::ValueFlag<std::string> argMovingFileName(allRequired, "moving", "The moving mesh (point-set) filename", {'m', "moving"});
  args::ValueFlag<std::vector<double>, args::DoubleVectorReader> argScale(allRequired, "scale", "The scale levels", {'s', "scale"});
  
  args::ValueFlag<std::string> argOutputFileName(parser, "output", "The output mesh (point-set) filename", {'o', "output"});
  args::ValueFlag<unsigned int> argNumberOfIterations(parser, "iterations", "The number of iterations", {'i', "iterations"});
  args::Flag argNormalize(parser, "normalize", "Normalization", {'n', "normalize"});
  args::Flag argTrace(parser, "trace", "Optimizer iterations tracing", {'T', "trace"});

  const std::string transformDescription =
    "The type of transform (That is number):\n"
    "  0 : Translation\n"
    "  1 : Versor3D\n"
    "  2 : Similarity\n"
    "  3 : ScaleSkewVersor3D\n";

  args::ValueFlag<size_t> argTypeOfTransform(parser, "transform", transformDescription, {'t', "transform"});
  
  const std::string metricDescription =
    "The type of metric (That is number):\n"
    "  0 : L2Rigid\n"
    "  1 : L2\n"
    "  2 : KC\n"
    "  3 : MLE\n";

  args::ValueFlag<size_t> argTypeOfMetric(parser, "metric", metricDescription, {'M', "metric"});

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help) {
    std::cout << parser;
    return EXIT_SUCCESS;
  }
  catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }
  catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }

  std::string fixedFileName = args::get(argFixedFileName);
  std::string movingFileName = args::get(argMovingFileName);
  std::vector<double> scale = args::get(argScale);
  
  unsigned int numberOfIterations = 100;

  if (argNumberOfIterations) {
    numberOfIterations = args::get(argNumberOfIterations);
  }

  bool normalize = argNormalize;
  bool trace = argTrace;

  size_t typeOfTransform = 0;

  if (argTypeOfTransform) {
    typeOfTransform = args::get(argTypeOfTransform);
  }

  size_t typeOfMetric = 0;

  if (argTypeOfMetric) {
    typeOfMetric = args::get(argTypeOfMetric);
  }

  std::cout << "options" << std::endl;
  std::cout << "number of iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  //
  MeshType::Pointer fixedMesh = MeshType::New();
  if (!readMesh<MeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << fixedFileName << std::endl;
  std::cout << "number of points " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  MeshType::Pointer movingMesh = MeshType::New();
  if (!readMesh<MeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << movingFileName << std::endl;
  std::cout << "number of points " << movingMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // read meshes
  if (!readMesh<MeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }
  PointSetType::Pointer movingPointSet = PointSetType::New();
  movingPointSet->SetPoints(movingMesh->GetPoints());

  if (!readMesh<MeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }
  PointSetType::Pointer fixedPointSet = PointSetType::New();
  fixedPointSet->SetPoints(fixedMesh->GetPoints());

  //--------------------------------------------------------------------
  // initialize scales
  size_t numberOfLevels = scale.size();
  itk::Array<double> fixedPointSetScale(scale.size());
  itk::Array<double> movingPointSetScale(scale.size());

  for (size_t n = 0; n < scale.size(); ++n) {
    fixedPointSetScale[n] = scale[n];
    movingPointSetScale[n] = scale[n];
  }

  // initialize transform
  typedef itk::PointSetPropertiesCalculator<PointSetType> PointSetPropertiesCalculatorType;
  PointSetPropertiesCalculatorType::Pointer fixedPointSetCalculator = PointSetPropertiesCalculatorType::New();
  fixedPointSetCalculator->SetPointSet(fixedPointSet);
  fixedPointSetCalculator->Compute();
  fixedPointSetCalculator->PrintReport(std::cout);

  PointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = PointSetPropertiesCalculatorType::New();
  movingPointSetCalculator->SetPointSet(movingPointSet);
  movingPointSetCalculator->Compute();
  movingPointSetCalculator->PrintReport(std::cout);

  typedef itk::Similarity3DTransform<double> InitialTransformType;
  InitialTransformType::Pointer fixedInitialTransform;
  InitialTransformType::Pointer movingInitialTransform;

  if (normalize) {
    fixedInitialTransform = InitialTransformType::New();
    fixedInitialTransform->SetIdentity();
    fixedInitialTransform->SetScale(1.0 / fixedPointSetCalculator->GetScale());
    fixedInitialTransform->SetCenter(fixedPointSetCalculator->GetCenter());
    InitialTransformType::TranslationType fixedTranslation;
    for (unsigned int n = 0; n < 3; ++n) {
      fixedTranslation[n] = -fixedPointSetCalculator->GetCenter()[n];
    }
    fixedInitialTransform->SetTranslation(fixedTranslation);

    movingInitialTransform = InitialTransformType::New();
    movingInitialTransform->SetIdentity();
    movingInitialTransform->SetScale(1.0 / movingPointSetCalculator->GetScale());
    movingInitialTransform->SetCenter(movingPointSetCalculator->GetCenter());
    InitialTransformType::TranslationType movingTranslation;
    for (size_t n = 0; n < 3; ++n) {
      movingTranslation[n] = -movingPointSetCalculator->GetCenter()[n];
    }
    movingInitialTransform->SetTranslation(movingTranslation);
  }
  else {
    for (size_t n = 0; n < scale.size(); ++n) {
      fixedPointSetScale[n] *= fixedPointSetCalculator->GetScale();
    }

    for (size_t n = 0; n < scale.size(); ++n) {
      movingPointSetScale[n] *= fixedPointSetCalculator->GetScale();
    }
  }

  // initialize transform
  typedef itk::InitializeTransform<double> InitializeTransformType;
  InitializeTransformType::Pointer initializeTransform = InitializeTransformType::New();
  if (!normalize) {
    initializeTransform->SetCenter(movingPointSetCalculator->GetCenter());
    initializeTransform->SetTranslation(fixedPointSetCalculator->GetCenter() - movingPointSetCalculator->GetCenter());
  }
  initializeTransform->SetTypeOfTransform(typeOfTransform);
  initializeTransform->Update();
  initializeTransform->PrintReport();
  TransformType::Pointer transform = initializeTransform->GetTransform();

  InitializeTransformType::ModeBoundsType modeBounds = initializeTransform->GetModeBounds();
  InitializeTransformType::BoundsType lowerBounds = initializeTransform->GetLowerBounds();
  InitializeTransformType::BoundsType upperBounds = initializeTransform->GetUpperBounds();

  //--------------------------------------------------------------------
  // initialize optimizer
  itk::LBFGSBOptimizer::Pointer optimizer = itk::LBFGSBOptimizer::New();
  optimizer->SetBoundSelection(modeBounds);
  optimizer->SetLowerBound(lowerBounds);
  optimizer->SetUpperBound(upperBounds);
  optimizer->SetMaximumNumberOfIterations(numberOfIterations);
  optimizer->MinimizeOn();
  if (trace) {
    typedef CommandIterationUpdate<itk::LBFGSBOptimizer> CommandLBFGSBOptimizerIterationUpdate;

    CommandLBFGSBOptimizerIterationUpdate::Pointer observer = CommandLBFGSBOptimizerIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);
  }

  //--------------------------------------------------------------------
  // metric
  typedef itk::InitializeMetric<PointSetType, PointSetType> InitializeMetricType;
  InitializeMetricType::Pointer initializeMetric = InitializeMetricType::New();
  initializeMetric->SetTypeOfMetric(typeOfMetric);
  try {
    initializeMetric->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializeMetric->PrintReport();
  InitializeMetricType::MetricType::Pointer metric = initializeMetric->GetMetric();

  //--------------------------------------------------------------------
  // perform registration
  typedef itk::GMMPointSetToPointSetRegistrationMethod<PointSetType, PointSetType> GMMPointSetToPointSetRegistrationMethodType;
  GMMPointSetToPointSetRegistrationMethodType::Pointer registration = GMMPointSetToPointSetRegistrationMethodType::New();
  registration->SetFixedPointSet(fixedPointSet);
  registration->SetFixedInitialTransform(fixedInitialTransform);
  registration->SetMovingPointSet(movingPointSet);
  registration->SetMovingInitialTransform(movingInitialTransform);
  registration->SetFixedPointSetScale(fixedPointSetScale);
  registration->SetMovingPointSetScale(movingPointSetScale);
  registration->SetOptimizer(optimizer);
  registration->SetMetric(metric);
  registration->SetTransform(transform);
  registration->SetNumberOfLevels(numberOfLevels);
  try {
    registration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }
  registration->Print(std::cout);

  std::cout << std::endl;
  std::cout << registration->GetMetric()->GetNameOfClass() << std::endl;
  std::cout << optimizer->GetStopConditionDescription() << std::endl;
  std::cout << "initial parameters " << registration->GetInitialTransformParameters() << std::endl;
  std::cout << "  final parameters " << optimizer->GetCurrentPosition() << std::endl;
  std::cout << "             value " << optimizer->GetValue() << std::endl;
  std::cout << std::endl;

  typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformMeshFilterType;
  TransformMeshFilterType::Pointer transformMesh = TransformMeshFilterType::New();
  transformMesh->SetInput(movingMesh);

  if (normalize) {
    typedef itk::CompositeTransform<double, PointSetType::PointDimension> CompositeTransformType;
    CompositeTransformType::Pointer transform = CompositeTransformType::New();
    transform->AddTransform(fixedInitialTransform->GetInverseTransform());
    transform->AddTransform(registration->GetTransform());
    transform->AddTransform(movingInitialTransform);
    transformMesh->SetTransform(transform);
  }
  else {
    transformMesh->SetTransform(transform);
  }

  try {
    transformMesh->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }

  if (argOutputFileName) {
    std::string fileName = args::get(argOutputFileName);
    std::cout << "write output mesh to the file " << fileName << std::endl;

    if (!writeMesh<MeshType>(transformMesh->GetOutput(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  PointSetType::Pointer outputPointSet = PointSetType::New();
  outputPointSet->SetPoints(transformMesh->GetOutput()->GetPoints());

  typedef itk::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
  PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(movingPointSet);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(outputPointSet);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  return EXIT_SUCCESS;
}
