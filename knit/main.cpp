
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string.h>
#include "CImg.h"
#include <omp.h>
#include <Windows.h>

using namespace cimg_library;

typedef int Chromosome;
typedef std::vector<Chromosome> DNA;
typedef std::pair<int, int> Point;
typedef CImg<float> Image;
typedef float DATA;

const int THREADS = 4;

double uniformRandom()
{
	return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}

double geometricRandom(double q) {
	double x = uniformRandom();
	return std::log(1 - x) / std::log(q);
}

double randomSign() {
	return rand() % 2 ? 1 : -1;
}

int intRandom(int start, int end) {
	return start + rand() % (end - start);
}
// return a normally distributed random number
double normalRandom()
{
	double u1 = uniformRandom();
	double u2 = uniformRandom();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

struct Triple {
	double min;
	double max;
	double mean;
};

int cross(const Point &O, const Point &A, const Point &B) {
	return (A.first - O.first) * (B.second - O.second) - (A.second - O.second) * (B.first - O.first);
}

bool operator <(const Point &a, const Point &b) {
	return a.first < b.first || (a.first == b.first && a.second < b.second);
}

bool isPointInPath(int x, int y, const std::vector<Point> &path) {
	int i = 0;
	int j = path.size() - 1;
	bool c = false;
	for (i = 0; i < path.size(); i++) {
		if ((path[i].second > y) != (path[j].second > y) && (x < (path[j].first - path[i].first) * (y - path[i].second) / (path[j].second - path[i].second) + path[i].first))
			c = !c;
		j = i;
	}
	return c;
}

bool isPointInCircle(int x, int y, int r) {
	int dx = r - x - 1;
	int dy = r - y - 1;
	return (dx * 2 + 1) * (dx * 2 + 1) + (dy * 2 + 1) * (dy * 2 + 1) <= (r * 2 + 1) * (r * 2 + 1);
}

void makeHull(const std::vector<Point> &nails, std::vector<Point> &bounds, int width, int height) {
	bounds.resize(height);
	int n = nails.size();


	/*std::vector<Point> points;
	points.assign(nails.begin(), nails.end());
	std::sort(points.begin(), points.end());
	std::vector<Point> H(2 * n);
	int k = 0;
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k - 2], H[k - 1], points[i]) < 0) k--;
		H[k++] = points[i];
	}

	// Build upper hull
	for (int i = n - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && cross(H[k - 2], H[k - 1], points[i]) < 0) k--;
		H[k++] = points[i];
	}

	H.resize(k - 1);*/

	for (int i = 0; i < height; i++) {
		bool wasInside = false;
		for (int j = 0; j < width; j++) {
			bool isInside = isPointInCircle(j, i, width / 2);
			if (isInside && !wasInside) {
				bounds[i].first = j;// == 0 ? 0 : j - 1;
				wasInside = true;
			}
			if (wasInside && !isInside) {
				bounds[i].second = j;
				wasInside = false;
			}
		}
		if (bounds[i].second == 0 && wasInside) {
			bounds[i].second = width - 1;
		}
	}
}

struct Slab {
	DATA *data;

	Slab() {
		data = new DATA[width * width];
		fill(0.0);
	}

	~Slab() {
		delete[] data;
	}

	void fill(DATA v = 0.0) {
		for (int i = 0; i < width * width; i++) {
			data[i] = v;
		}
	}

	void fromImage(const Image &image) {
		for (int i = 0; i < image.width(); i++) {
			for (int j = 0; j < image.height(); j++) {
				if (i >= bounds[j].first && i <= bounds[j].second) {
					data[i + j * width] = image(i, j, 0, 0);
				}
			}
		}
	}

	void toImage(Image &image) {
		image.fill(0.0);
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				image(j, i, 0, 0) = data[j + i * width];
			}
		}
	}

	void normalize() {
		DATA mean = 0.0;
		DATA dev = 0.0;
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				mean += data[j + i * width];
			}
		}
		mean /= size;

		DATA d;
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				d = data[j + i * width] - mean;
				dev += d * d;
			}
		}
		dev = sqrt(dev);

		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				data[j + i * width] = (data[j + i * width] - mean) / dev;
			}
		}
	}

	double covariate(const Slab &kernel) {
		DATA result = 0.0;
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				result += data[j + i * width] * kernel.data[j + i * width];
			}
		}
		return result;
	}

	void addPixel(int x, int y, DATA opacity) {
		data[y * width + x] += opacity;
	}

	void clamp() {
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				data[j + i * width] = std::min(data[j + i * width], (DATA)1.0);
			}
		}
	}

	void pow(DATA power) {
		for (int i = 0; i < height; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				data[j + i * width] = std::pow(data[j + i * width], power);
			}
		}
	}

	void scan() {
		int j = height / 2;
		for (int i = 0; i < width; i++) {
			printf("%i:\t%.7lf\n", i, data[j * width + i]);
		}
	}
	
	static int width;
	static int height;
	static int size;
	static std::vector<Point> bounds;
	static void prepare(int imageWidth, int imageHeight, const std::vector<Point> nails) {
		width = imageWidth;
		height = imageHeight;
		bounds.resize(height);
		makeHull(nails, bounds, width, height);
		int l = 0;
		for (int i = 0; i < width; i++) {
			l += bounds[i].second - bounds[i].first;
		}
		size = l;
	}
};

int Slab::width;
int Slab::height;
int Slab::size;
std::vector<Point> Slab::bounds;

struct Genome {
	DNA dna;
	double fitness;
	double weight;

	void init(int size, int nails) {
		dna.resize(size);
		int minOffset = nails * 1 / 5;
		int maxOffset = nails * 4 / 5;
		dna[0] = rand() % nails;
		for (int i = 1; i < size; i++) {
			dna[i] = rand() % nails;    //(dna[i - 1] + intRandom(minOffset, maxOffset)) % nails;
		}
		fitness = 0.0;
	}

	void mutate(double count, double spread, int nails) {
		for (int i = 0; i < dna.size(); i++) {
			double x = uniformRandom();
			if (x < count) {
				x = normalRandom();
				int offset = floor(x * spread * nails + 0.5);
				while (offset < 0) offset += nails; // to eliminate negative offsets and %
				dna[i] = (dna[i] + offset) % nails;
			}
		}
	}

	void burst(float probability, int nails) {
		for (int i = 0; i < dna.size(); i++) {
			if (uniformRandom() < probability) {
				int offset = floor(uniformRandom() * nails);
				dna[i] = offset;
			}
		}
	}

	void draw(Slab &canvas, const std::vector<Point> &nails, const float threadOpacity) {
		canvas.fill(0.0);
		DATA paint = threadOpacity;
		int x0, x1, y0, y1, dx, dy, sx, sy, err;
		for (int i = 0; i < dna.size() - 1; i++) {
			x0 = nails[dna[i]].first;
			x1 = nails[dna[i + 1]].first;
			y0 = nails[dna[i]].second;
			y1 = nails[dna[i + 1]].second;
			dx = cimg::abs(x1 - x0);
			dy = cimg::abs(y1 - y0);
			sx = (x0 < x1) ? 1 : -1;
			sy = (y0 < y1) ? 1 : -1;
			err = dx - dy;
			while (true) {
				canvas.addPixel(x0, y0, paint);			

				if ((x0 == x1) && (y0 == y1)) break;
				int e2 = 2 * err;
				if (e2 >-dy) { err -= dy; x0 += sx; }
				if (e2 < dx) { err += dx; y0 += sy; }
			}
		}
		canvas.clamp();
	}

	// ORDER BY fitness DESC
	bool operator <(const Genome &g) {
		return fitness > g.fitness;
	}
};

typedef std::vector<Genome> Genomes;

struct Population {
	Genomes population;
	double weightSum;

	void init(int populationSize, int genomeSize, int nails) {
		population.resize(populationSize);
		for (int i = 0; i < population.size(); i++) {
			population[i].init(genomeSize, nails);
		}
	}

	double calculateFitness(const Slab &kernel, Slab* canvases, const std::vector<Point> &nails, const float threadOpacity) {
		double minFitness = 10e+7;
		double maxFitness = -10e+7;
		double average = 0.0;
		int size = population.size();

#pragma omp parallel for
		for (int i = 0; i < size; i++) {
			int this_thread = omp_get_thread_num();
			//printf("%i\n", i);
			population[i].draw(canvases[this_thread], nails, threadOpacity);
			//canvases[this_thread].pow(0.75);
			canvases[this_thread].normalize();
			double fitness = canvases[this_thread].covariate(kernel);
			population[i].fitness = fitness;
			if (fitness > maxFitness) maxFitness = fitness;
			if (fitness < minFitness) minFitness = fitness;
			average += fitness;
		}
		std::sort(population.begin(), population.end());

		weightSum = 0.0;
		for (int i = 0; i < population.size(); i++) {
			double weight = (population[i].fitness - minFitness) / (maxFitness - minFitness);
			population[i].weight = weightSum;
			weightSum += weight;
		}
		
		return average / population.size();
	}

	Genome* getParent() {
		double random = uniformRandom() * weightSum / 2; // only first (successful) part, because fuck you
		int i;
		for (i = population.size() - 1; i > 0; i--) {
			if (random > population[i].weight)
				break;
		}
		//printf("Parent: %i (%lf)\n", i, population[i].fitness);
		return &(population[i]);
	}

	Genome* getParentOrdinal() {
		double random = pow(uniformRandom(), 2) * population.size() / 2;
		int i = floor(random);
		//printf("%i ", i);
		return &(population[i]);
	}

	void burst(float probability, int nails) {
		for (int i = 0; i < population.size(); i++) {
			population[i].burst(probability, nails);
		}
	}
};

double MAX_CHUNK, MIN_CHUNK;
int getCrossoverJumpSize(int size) {
	return floor(size * ( uniformRandom() * (MAX_CHUNK - MIN_CHUNK) + MIN_CHUNK) );
}

Genome crossover(Genome *mother, Genome *father) {
	Genome child;
	int size = mother->dna.size();

	child.dna.resize(size);

	const Genome *donor;
	int i = 0;
	int j = 0;
	int jump = getCrossoverJumpSize(size);
	donor = uniformRandom() < 0.5 ? mother : father;
	while (i < size) {
		child.dna[i] = donor->dna[i];
		i++;
		j++;
		if (j == jump) {
			j = 0;
			jump = getCrossoverJumpSize(size);
			donor = donor == father ? mother : father;
		}
	}
	return child;
}

void generateNewPopulation(Population &base, Population &result, int nails, int elite, double mutationCount, double mutationSpread) {
	for (int i = 0; i < elite && i < result.population.size(); i++) {
		result.population[i] = base.population[i];
	}
	int size = result.population.size();
#pragma omp parallel for
	for (int i = elite; i < size; i++) {
		// select 2 parents
		Genome* mother = base.getParentOrdinal();
		Genome* father = base.getParentOrdinal();

		// crossover
		result.population[i] = crossover(mother, father);

		// mutate
		result.population[i].mutate(mutationCount, mutationSpread, nails);
	}
}

Image loadImage(char *filename, int imageSize) {
	Image result(imageSize, imageSize, 1, 1, 0.0);
	CImg<unsigned char> image(filename);

	image.resize(imageSize, imageSize);

	for (int i = 0; i < imageSize; i++) {
		for (int j = 0; j < imageSize; j++) {
			double c = image(i, j, 0, 0);
			c = 1.0 - c / 255.0;
			result(i, j, 0, 0) = c;
		}
	}

	return result;
}

const int PROCESS_MODE = 0;
const int DRAW_MODE = 1;

void drawComplete(FILE *data, int nails, int imageSize, double opacity, const std::vector<Point> &nailPoints) {
	int size;
	double threadDiameter = 0.0, ringDiameter = 0.0;
	fscanf_s(data, "%i %lf %lf", &size, &threadDiameter, &ringDiameter);

	opacity = imageSize * threadDiameter / ringDiameter;

	printf("Thread opacity: %lf\n", opacity);

	Genome gen;
	gen.init(size, nails);

	for (int i = 0; i < size - 1; i++) {
		int t;
		fscanf_s(data, "%i", &t);
		gen.dna[i] = t;
	}

	Slab canvas;
	gen.draw(canvas, nailPoints, opacity);

	Image output(imageSize, imageSize, 1, 1, 0.0f);
	canvas.toImage(output);

	output.normalize(0, 255);
	output *= -1;
	output += 255;

	output.display();
}

BOOL ctrlHandler(DWORD fdwCtrlType) {
	return TRUE;
}

void run(char *datafile) {
	int nails = 150;
	int populationSize = 100;
	int elite = 0;
	int genomeSize = nails * 4;
	int imageSize = 100;
	double minChunk = 0.5;
	double maxChunk = 0.5;
	int iterations = 50;
	int stableIterations = 5;
	int iterationsToBurst = 100;
	float threadOpacity = 1.0 / 10;
	double mutationSpread = 1.0;
	double mutationCount = 1.0;
	double burstProbability = 0.005;
	double threshold = 1e-7;
	double ringDiameter = 0;
	double threadDiameter = 0;
	int saveEvery = -1;
	int saveByFitness = -1;
	int historyStep = 1;
	char filename[256] = "test-i.png";

	int mode = PROCESS_MODE;
	FILE *config = NULL;
	if (datafile != NULL) {
		mode = DRAW_MODE;
		fopen_s(&config, datafile, "r");
	} else {
		fopen_s(&config, "config.txt", "r");
	}
	
	
	char line[256];
	while (fscanf_s(config, "%s", line, 256) != EOF) {
		if (strcmp(line, "file") == 0) {
			//strcpy_s(filename, line);
			fscanf_s(config, "%s", filename, 256);
		} else
		if (strcmp(line, "nails") == 0) {
			fscanf_s(config, "%i", &nails);
		} else
		if (strcmp(line, "populationsize") == 0) {
			fscanf_s(config, "%i", &populationSize);
		} else
		if (strcmp(line, "elite") == 0) {
			fscanf_s(config, "%i", &elite);
		} else
		if (strcmp(line, "genomesize") == 0) {
			fscanf_s(config, "%i", &genomeSize);
		} else
		if (strcmp(line, "minchunk") == 0) {
			fscanf_s(config, "%lf", &minChunk);
			MIN_CHUNK = minChunk;
		} else
		if (strcmp(line, "maxChunk") == 0) {
			fscanf_s(config, "%lf", &maxChunk);
			MAX_CHUNK = maxChunk;
		} else
		if (strcmp(line, "imagesize") == 0) {
			fscanf_s(config, "%i", &imageSize);
		} else
		if (strcmp(line, "iterations") == 0) {
			fscanf_s(config, "%i", &iterations);
		}
		if (strcmp(line, "mutation") == 0) {
			fscanf_s(config, "%lf", &mutationSpread);
		} else
		if (strcmp(line, "mutationcount") == 0) {
			fscanf_s(config, "%lf", &mutationCount);
		}
		else
		if (strcmp(line, "stableiterations") == 0) {
			fscanf_s(config, "%i", &stableIterations);
		} else
		if (strcmp(line, "goodgenerations") == 0) {
			fscanf_s(config, "%i", &iterationsToBurst);
		} else
		if (strcmp(line, "threshold") == 0) {
			fscanf_s(config, "%lf", &threshold);
		} else
		if (strcmp(line, "opacity") == 0) {
			fscanf_s(config, "%lf", &threadOpacity);
		} else
		if (strcmp(line, "burst") == 0) {
			fscanf_s(config, "%lf", &burstProbability);
		}
		if (strcmp(line, "ring_diameter") == 0) {
			fscanf_s(config, "%lf", &ringDiameter);
		} else
		if (strcmp(line, "save_every") == 0) {
			fscanf_s(config, "%i", &saveEvery);
		} else
		if (strcmp(line, "save_by_fitness") == 0) {
			fscanf_s(config, "%i", &saveByFitness);
		}
		else
		if (strcmp(line, "history_step") == 0) {
			fscanf_s(config, "%i", &historyStep);
		} else
		if (strcmp(line, "thread_diameter") == 0) {
			fscanf_s(config, "%lf", &threadDiameter);
		} else 
		if (strcmp(line, "---") == 0) {
			break;
		} else {
			
		}
	}

	if (mode == DRAW_MODE) {
		imageSize = 500;
	}

	if (abs(ringDiameter) > 1e-7 && abs(threadDiameter) > 1e-7) {
		threadOpacity = imageSize * threadDiameter / ringDiameter;
	}

	std::vector<Point> nailPoints;
	nailPoints.resize(nails);
	for (int i = 0; i < nails; i++) {
		double a = 3.14159265359 * 2 * i / nails;
		double s = std::sin(a);
		double c = std::cos(a);
		double r = imageSize / 2 - 1;
		nailPoints[i].first = floor(imageSize / 2 + c * r);
		nailPoints[i].second = floor(imageSize / 2 + s * r);
	}

	Slab::prepare(imageSize, imageSize, nailPoints);

	if (mode == DRAW_MODE) {
		drawComplete(config, nails, imageSize, threadOpacity, nailPoints);
		fclose(config);
		return;
	}

	fclose(config);

	//SetConsoleCtrlHandler()

	Population populations[2];
	populations[0].init(populationSize, genomeSize, nails);
	populations[1].init(populationSize, genomeSize, nails);

	Image kernelImage = loadImage(filename, imageSize);
	Image output(imageSize, imageSize, 1, 1, 0.0f);

	Slab kernel;
	kernel.fromImage(kernelImage);
	kernel.normalize();

	/*kernel.toImage(kernelImage);
	kernelImage.display();*/

	Slab canvas;

	omp_set_dynamic(0);
	omp_set_num_threads(THREADS);
	Slab *canvases = new Slab[THREADS];

	double fitness = 0.0;
	double prevFitness = 1.0;
	double avgTime = 0.0;
	const int timeDepth = 50;
	std::vector<double> timeSeries;
	timeSeries.resize(timeDepth);

	double best = -1.0;
	double prevBest = -1.0;
	int stableGeneration = 0;
	int goodGeneration = 0;
	int generation = 0;
	bool enough = false;
	std::vector<Triple> history;
	history.reserve(iterations);
	
	int intFitness = -1;

	while (generation < iterations) {
		int current = generation % 2;
		int next = (generation + 1) % 2;

		time_t start = clock();

		prevFitness = fitness;
		fitness = populations[current].calculateFitness(kernel, canvases, nailPoints, threadOpacity);

		prevBest = std::max(prevBest, best);
		best = populations[current].population.front().fitness;

		Triple h;
		h.max = populations[current].population.front().fitness;
		h.min = populations[current].population.back().fitness;
		h.mean = fitness;
		history.push_back(h);

		int currentIntFitness = (int)floor(h.max * 100);

		if ((saveEvery > 0 && generation % saveEvery == 0) || (saveByFitness > 0 && currentIntFitness > intFitness && currentIntFitness % saveByFitness == 0)) {
			//printf("SAVING\n");
			Genome *win = &(populations[generation % 2].population.front());
			win->draw(canvas, nailPoints, threadOpacity);
			canvas.toImage(output);
			char fname[100];
			sprintf_s(fname, 100, "temp\\generation-%i.bmp", generation);
			output.normalize(0, 255);
			output *= -1;
			output += 255;
			output.save(fname);
		}

		intFitness = currentIntFitness;

		generateNewPopulation(populations[current], populations[next], nails, elite, mutationCount, mutationSpread);

		start = clock() - start;
		double seconds = ((double)start / CLOCKS_PER_SEC);

		if (generation == 0) {
			for (int i = 0; i < timeSeries.size(); i++) {
				timeSeries[i] = seconds;
			}
		}

		timeSeries[generation % timeDepth] = seconds;
		avgTime = 0;
		for (int i = 0; i < timeSeries.size(); i++) {
			avgTime += timeSeries[i];
		}
		avgTime /= timeSeries.size();
		printf("G: %i (%.2lf%%), F: %.7lf, T: %.5lfs (E: %.2lfm)\r", generation, (double)generation * 100 / iterations, h.max, avgTime, (iterations - generation) * avgTime / 60);

		if (best < prevBest + 1e-10) {
			goodGeneration++;
		} else {
			goodGeneration = 0;
		}

		if (goodGeneration > iterationsToBurst) {
			goodGeneration = 0;
			populations[next].burst(burstProbability, nails);
			//printf("BUSRT!\n");
		}

		if (abs(1.0 - fitness / prevFitness) < threshold) {
			stableGeneration++;
		} else {
			stableGeneration = 0;
		}
		if (stableGeneration >= stableIterations) {
			enough = true;
		}
		generation++;
	}
	fitness = populations[generation % 2].calculateFitness(kernel, canvases, nailPoints, threadOpacity);
	printf("Generation: %i,\tFitness: %lf\n", generation, fitness);

	Genome *win = &(populations[generation % 2].population.front());
	win->draw(canvas, nailPoints, threadOpacity);

	canvas.toImage(output);

	output.normalize(0, 255);
	output *= -1;
	output += 255;

	double length = 0.0;
	for (int i = 0; i < win->dna.size() - 1; i++) {
		Point p1 = nailPoints[win->dna[i]];
		Point p2 = nailPoints[win->dna[i+1]];
		int dx = p1.first - p2.first;
		int dy = p1.second - p2.second;
		length += sqrt(dx * dx + dy * dy);
	}
	length *= ringDiameter / imageSize;
	printf("Length: %.2lf m\n", length / 1000.0);

	output.save("result.bmp");

	output.display("Winner");

	FILE *outfile;
	fopen_s(&outfile, "history.csv", "w");
	for (int i = 0; i < history.size(); i++) {
		if(i % historyStep == 0)
			fprintf_s(outfile, "%i\t%.7lf\t%.7lf\t%.7lf\n", i, history[i].min, history[i].mean, history[i].max);
	}
	fclose(outfile);

	fopen_s(&outfile, "result.txt", "w");
	fprintf_s(outfile, "file %s\nnails %i\nimagesize %i\npopulationsize %i\nelite %i\ngenomesize %i\nmutationcount %lf\nmutation %lf\nminchunk %lf\nmaxchunk %lf\niterations %i\nstableiterations %i\ngoodgenerations %i\nburst %lf\nthreshold %lf\nring_diameter %lf\nthread_diameter %lf\nsave_every %i\nhistory_step %i",
		filename, nails, imageSize, populationSize, elite, genomeSize, mutationCount, mutationSpread, minChunk, maxChunk, iterations, stableIterations, iterationsToBurst, burstProbability, threshold, ringDiameter, threadDiameter, saveEvery, historyStep);
	fprintf_s(outfile, "\n---\n");
	fprintf_s(outfile, "%i %lf %lf ", genomeSize, threadDiameter, ringDiameter);
	for (int i = 0; i < win->dna.size() - 1; i++) {
		fprintf_s(outfile, " %i", win->dna[i]);
	}
	fprintf_s(outfile, "\n%lf meters", length / 1000.0);
	fclose(outfile);

	delete[] canvases;
	//drawHistory(history);
}

int main(int argc, char **argv) {
	srand(clock());

	/*int vals[100];
	for (int i = 0; i < 100; i++) vals[i] = 0;
	int n = 1000;
	for (int i = 0; i < n; i++) {
		int k = floor(normalRandom() * 200 * 0.05);
		printf("%i\n", k);
	}*/

	if (argc == 1)
		run(NULL);
	else
		run(argv[1]);

	char c;
}