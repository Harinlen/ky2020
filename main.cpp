#define _CRT_SECURE_NO_WARNINGS

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <queue>

#include <omp.h>

using namespace std;
using std::priority_queue;

// Macros.
#define            N                 (30)
#define            ID_VERTEX         (0)
#define            ID_VALUE          (1)
const float        EPSINON          = 0.00001f;
// Constant values.
static const float CONNECT          = 1.0f;
static const int   EDGE_PARAMS      = 35;
static const float I                = 1073741824.0f;
// File paths.
static const char *PATH_INPUT_EDGE  = "edges";
static const char *PATH_NETWORK     = "network";
static const char *PATH_OUTPUT      = "output";
// Public shared variables.
static float **    results          = nullptr;
static float **    edge_value       = nullptr;
static int         edges            = 0;
static int         vertex[N][N][2];
static int         cost_mapping[N]  = {
    -1, // 0
    -1, // 1
    -1, // 2
    0,  // 3
    -1, // 4
    34, // 5
    3,  // 6
    2,  // 7
    14, // 8
    12, // 9
    -1, // 10
    29, // 11
    15, // 12
    -1, // 13
    10, // 14
    26, // 15
    24, // 16
    27, // 17
    23, // 18
    -1, // 19
    -1, // 20
    -1, // 21
    -1, // 22
    -1, // 23
    -1, // 24
    31, // 25
    19, // 26
    -1, // 27
    -1, // 28
    1   // 29
};
//Priority queue structure.
typedef struct QueueNode
{
    int vid;
    float value;
    QueueNode(int mvid=0, float mvalue=0) :
        vid(mvid),
        value(mvalue)
    {}

    bool operator < (const QueueNode& rhs) const {
      return value > rhs.value;
    }
} QueueNode;

bool is_close(const float &f1, const float &f2)
{
    const float diff = f1 - f2;
    return (diff >= -EPSINON) && (diff <= EPSINON);
}

FILE *open_file_read(char const* filePath, const char *errorMessage)
{
    FILE *file = fopen(filePath, "r");
    if(file == nullptr)
    {
        printf("%s\n", errorMessage);
        exit(-1);
    }
    return file;
}

FILE *open_file_write(char const* filePath, const char *errorMessage)
{
    FILE *file = fopen(filePath, "w");
    if(file == nullptr)
    {
        printf("%s\n", errorMessage);
        exit(-1);
    }
    return file;
}

void load_edges(const char *edge_path = nullptr)
{
    // Read the edge value file.
    FILE *edge_file = nullptr;
    if(edge_path)
    {
        edge_file = open_file_read(edge_path,
                                   "cannot open edge value file");
    }
    else
    {
        edge_file = open_file_read(PATH_INPUT_EDGE,
                                   "cannot open edge value file");
    }
    // Read the line of the edges.
    fscanf(edge_file, "%d", &edges);
    // Create edge buffers.
    float *edge_buffer = static_cast<float *>(malloc(sizeof(float) * static_cast<size_t>(edges * EDGE_PARAMS)));
    // Create edges array.
    edge_value = static_cast<float **>(malloc(sizeof(float *) * static_cast<size_t>(edges)));
    for(int i=0; i<edges; ++i)
    {
        // Create edge buffer.
        edge_value[i] = edge_buffer + i * EDGE_PARAMS;
        // Read edge data.
        for(int j=0; j<EDGE_PARAMS; ++j)
        {
            fscanf(edge_file, "%f", &edge_value[i][j]);
        }
    }
    // Close file.
    fclose(edge_file);
}

void load_cost(const char *network_path = nullptr)
{
    /* Index equals list:
     * 0, Positive: the index of the edge list.
     * -1: Constant value 0.
     * -2: Constant value I.
     * -3: Constant value CONNECT.
     */
    // Read the cost file.
    FILE *network_file = open_file_read(network_path ? network_path : PATH_NETWORK,
                                        "cannot open network input file");
    // Create the cost network.
    for(int i=0; i<N; ++i)
    {
        // Clear the counter.
        vertex[i][0][ID_VERTEX] = 0;
    }
    // Read the network file until it reaches the end of the file.
    while(!feof(network_file))
    {
        int i, j;
        // Read one line.
        if(fscanf(network_file,"%d%d", &i, &j) == 2)
        {
            // Increase the counter.
            ++vertex[i-1][0][ID_VERTEX];
            // Save the connected index.
            vertex[i-1][vertex[i-1][0][ID_VERTEX]][ID_VERTEX] = j-1;
            // Save the connected value.
            vertex[i-1][vertex[i-1][0][ID_VERTEX]][ID_VALUE] =
                    (cost_mapping[i] == -1) ? -3 : cost_mapping[i];
        }
    }
    // Sync all the counter index.
    for(int i=0; i<N; ++i)
    {
        // Clear the counter.
        vertex[i][0][ID_VALUE] = vertex[i][0][ID_VERTEX];
    }
    fclose(network_file);
}

void dijkstra(float *multi, float *results)
{
    // Prepare the variables.
    bool final[N];
    float dist[N];
    QueueNode current;
    // Prepare the NE and NF.
    float NE = 0.0f, NF = 0.0f;
    // Initial the distance and final state.
    for(int m=0; m<N; ++m)
    {
        // Execute Dijkstra for all vertex.
        for(int i=0; i<N; ++i)
        {
            final[i] = false;
            dist[i] = (i == m) ? 0.0f : I;
        }
        // Prepare the priority queue.
        priority_queue<QueueNode> queue;
        current.vid = m;
        current.value = 0.0f;
        // Put the current node to queue.
        queue.push(current);
        // Loop until all the queue is executed.
        while(!queue.empty())
        {
            // Pop the top from the queue.
            current = queue.top();
            queue.pop();
            // Check final, although I think it is unnecessary.
            if(final[current.vid])
            {
                // Do not need to search again.
                continue;
            }
            // Mark as final.
            final[current.vid] = true;
            // Loop for checking all the edges.
            for(int i=1; i<=vertex[current.vid][0][ID_VERTEX]; ++i)
            {
                int to_id = vertex[current.vid][i][ID_VERTEX];
                // Check final of the vid point.
                if(!final[to_id])
                {
                    int edge_id = vertex[current.vid][i][ID_VALUE];
                    float edge = (edge_id == -3) ? CONNECT : multi[edge_id],
                          trial = edge + dist[current.vid];
                    if(dist[to_id] > trial)
                    {
                        // Update the distance.
                        dist[to_id] = trial;
                        // Add the edge to queue.
                        queue.push(QueueNode(to_id, dist[to_id]));
                    }
                }
            }
        }
        // Calculate NE and NF.
        // NE
        for(int i=0; i<N; ++i)
        {
            if(!is_close(dist[i], 0) && !is_close(dist[i], I))
            {
                NE += 1.0f / dist[i];
            }
        }
        // NF.
        if(m!=29)
        {
            NF += 1.0f / dist[29];
        }
    }
    // Save the result.
    results[0] = NE;
    results[1] = NF;
}

int main(int argc, char *argv[])
{
    double time_start, time_end, time_cost;
    time_start=clock();
    // --- Start ---
    // Load edge values.
    load_edges((argc > 1) ? argv[1] : nullptr);
    // Load cost.
    load_cost((argc > 2) ? argv[2] : nullptr);
    // Create result buffer.
    float *result_buffer = static_cast<float *>(malloc(sizeof(float) * static_cast<size_t>(edges) * 2));
    results = static_cast<float **>(malloc(sizeof(float *) * static_cast<size_t>(edges)));
    for(int i=0; i<edges; ++i)
    {
        results[i] = result_buffer + i * 2;
    }
    // Loop and calculate for each edge.
    #pragma omp parallel for
    for(int i=0; i<edges; ++i)
    {
        dijkstra(edge_value[i], results[i]);
    }
    // Open the target file.
    FILE *output_file = open_file_write((argc > 3) ? argv[3] : PATH_OUTPUT,
                                        "cannot open output file");
    for(int i=0; i<edges; ++i)
    {
        fprintf(output_file, "%8.6f\t%8.6f\n", results[i][0], results[i][1]);
    }
    fclose(output_file);
    // --- Finish ---
    time_end=clock();
    time_cost=time_end-time_start;
    printf("%f\n",time_cost);
    return 0;
}
