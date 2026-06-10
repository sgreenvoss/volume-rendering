#include <vtkImageData.h>
#include <vtkStaticCellLocator.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPNGWriter.h>
#include <Kokkos_Core.hpp>
#include <vtkObject.h>
#include <vtkGenericCell.h>
#include <vtkSMPTools.h>
#include <omp.h>
#include <chrono>
#include <fstream>

#include "DataStructs.h"

// parallel x, threshold 0.99: 402 seconds
// parallel x, threshold 0.95: 643 seconds?
// parallel x, threshold 0.95: (attempt 2) 374 seconds?
// parallel x, threshold 0.90: 374 seconds?
// more more parallel: 852 :( 791 669

// properly enabling omp, threshold 0.95: 200s
// parallel ray sampling too, threshold 0.95: 314
// parallel ray sampling too, removing shared vars smh, threshold 0.95: 319

using namespace std;

const double PI = 3.1415926;

double colorEq(double front_o, double back_o, double front, double back) {
	return front + (1 - front_o) * back_o * back;
	// 0.0000092 + (.002355) *  0.304082 * 0
}

double map225to1(unsigned char val) {
	int v = val;
	double slope = 1.0 / 255.0;
	/*cout << "v: " << v << endl;
	cout << "double v: " << (double)v << endl;
	cout << "mapped color " << v << " to " << slope * (double) v << endl;*/
	return slope * (double) v;
}

unsigned char map1to255(double val) {
	double slope = 255.0;
	return (unsigned char)(slope * val);
}

void sampleAlongRay(Vector3<double> ray, const int n_samples, 
					Camera cam, int step_size, double* samples,
					vtkRectilinearGrid* data, vtkStaticCellLocator* locator,
					vtkDataArray* scalars) {

	vtkGenericCell* genericCell = vtkGenericCell::New();

	const double ox = cam.position.x + ray.x * cam.near;
	const double oy = cam.position.y + ray.y * cam.near;
	const double oz = cam.position.z + ray.z * cam.near;
	const double dx = ray.x * step_size;
	const double dy = ray.y * step_size;
	const double dz = ray.z * step_size;

	for (int i = 0; i < n_samples; i++) {
		double pcoords[3];
		double weights[8];
		double point[3];
		point[0] = ox + i * dx;
		point[1] = oy + i * dy;
		point[2] = oz + i * dz;
		
		vtkIdType cellId = locator->FindCell(point, 1e-10, genericCell, pcoords, weights);

		if (cellId < 0) {
			samples[i] = 0.0;
			continue;
		}

		double interp_val = 0.0;
		for (int j = 0; j < 8; j++) {
			interp_val += weights[j] * scalars->GetTuple1(genericCell->GetPointId(j));
		}
		samples[i] = interp_val;
	}
	genericCell->Delete();
}

//struct ProcessRayFunctor {
//	double* _image_buffer;
//
//	ProcessRayFunctor(double* )
//};


int main(int argc, char* argv[])
{
	vtkObject::GlobalWarningDisplayOff();
	//vtkSMPTools::SetBackend("OpenMP");
	cout << _OPENMP << endl;
	const int IMG_SIZE = 1000;
	const int REF_SAMPLE = 500;
	const int SAMPLES_PER_RAY = 500;
	double op_ratio = (double)REF_SAMPLE / (double)SAMPLES_PER_RAY;
	const char* FILE_NAME = "astro512_ascii.vtk";

	auto image = vtkImageData::New();
	image->SetDimensions(IMG_SIZE, IMG_SIZE, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	unsigned char* image_buffer = (unsigned char*) image->GetScalarPointer();

	ifstream inputFile(FILE_NAME);
	Camera cam = SetupCamera();
	TransferFunction tf = SetupTransferFunction();

	// vector initialization for finding rays
	Vector3<double> look = cam.focus - cam.position;
	Vector3<double> look_norm = look / look.Magnitude();

	Vector3<double> lxu = look.Cross(cam.up);
	Vector3<double> ru = lxu / lxu.Magnitude();
	Vector3<double> lxru = look.Cross(ru);
	Vector3<double> rv = lxru / lxru.Magnitude();
	double rads = cam.angle * (PI / 180.0);
	double k = (2.0 * tan(rads / 2.0)) / (double)IMG_SIZE;
	Vector3<double> dx = k * ru;
	Vector3<double> dy = k * rv;

	// getting extents to calculate step size
	vtkSmartPointer<vtkRectilinearGridReader> reader =
		vtkSmartPointer<vtkRectilinearGridReader>::New();
	reader->SetFileName(FILE_NAME); reader->Update();
	vtkRectilinearGrid* data = reader->GetOutput();

	double step_size = (cam.far - cam.near) / (double)(SAMPLES_PER_RAY - 1);
	cout << "step size: " << step_size << endl;

	cout << "ru: " << ru << endl;
	cout << "rv: " << rv << endl;
	cout << "dx " << dx << endl;
	cout << "dy: " << dy << endl;

	vtkPointData* pt_data = data->GetPointData();
	string arr_name = pt_data->GetArrayName(0);
	vtkDataArray* scalars = pt_data->GetArray(arr_name.c_str());

	auto start = chrono::steady_clock::now();
	cout << "started the for loop" << endl;

	vtkStaticCellLocator* locator = vtkStaticCellLocator::New();
	locator->SetDataSet(data);
	locator->BuildLocator();

#pragma omp parallel
{
	unsigned char* local_colors = tf.colors;
	double* local_opacities = tf.opacities;
	const int local_num_bins = tf.numBins;
	const double local_min = tf.min;
	const double local_max = tf.max;

	#pragma omp for collapse(2)
	for (int y = 0; y < IMG_SIZE; y++) {
			for (int x = 0; x < IMG_SIZE; x++) {
				Vector3<double> ray = look_norm + ((2.0 * (double)x + 1.0 - (double)IMG_SIZE) / 2.0) * dx
					+ ((2.0 * (double)y + 1.0 - (double)IMG_SIZE) / 2.0) * dy;

				double sample[SAMPLES_PER_RAY];

				sampleAlongRay(ray, SAMPLES_PER_RAY, cam, step_size, sample, data, locator, scalars);

				// calculate color
				double sampleRGB[3];
				double final_opacity;

				for (int i = 0; i < SAMPLES_PER_RAY; i++) {
					unsigned char RGB[3];
					double opacity;
					int bin = TransferFunction::GetBin(sample[i], local_max, local_min, local_num_bins);
					TransferFunction::ApplyTransferFunction(sample[i], RGB, opacity, local_colors, bin, local_opacities);

					if (i >= 1) {
						opacity = 1 - pow(1 - opacity, op_ratio);
						sampleRGB[0] = colorEq(final_opacity, opacity, sampleRGB[0], map225to1(RGB[0]));
						sampleRGB[1] = colorEq(final_opacity, opacity, sampleRGB[1], map225to1(RGB[1]));
						sampleRGB[2] = colorEq(final_opacity, opacity, sampleRGB[2], map225to1(RGB[2]));
						final_opacity = final_opacity + (1 - final_opacity) * opacity;
					}
					else {
						sampleRGB[0] = map225to1(RGB[0]);
						sampleRGB[1] = map225to1(RGB[1]);
						sampleRGB[2] = map225to1(RGB[2]);
						final_opacity = opacity;
					}

					if (final_opacity >= 0.99) {
						break;
					}
				}

				int img_index = (y * IMG_SIZE + x) * 3;
				image_buffer[img_index] = map1to255(sampleRGB[0]);
				image_buffer[img_index + 1] = map1to255(sampleRGB[1]);
				image_buffer[img_index + 2] = map1to255(sampleRGB[2]);
			}
	}
}
	auto end = std::chrono::steady_clock::now();
	auto duration = end - start; 

	auto seconds = chrono::duration_cast<chrono::seconds>(duration).count();
	cout << "execution time: " << seconds << " seconds" << endl;

	auto writer = vtkPNGWriter::New();
	writer->SetFileName("99astro.png");
	writer->SetInputData(image);
	writer->Write();

}


