using System;
using System.IO;
using GMath;
using static GMath.Gfx;
using static GMath.GRandom;
using static GMath.GTools;

namespace SphereScattering
{
    class Program
    {
        /// <summary>
        /// Gets the variables for a final state of a path in a spherical medium
        /// </summary>
        public class PathSummary
        {
            /// <summary>
            /// Number of scatters. Allways >= 1 because the center produces the first scatter
            /// </summary>
            public int N;
            /// <summary>
            /// Final position in the surface of the unitary sphere
            /// </summary>
            public float3 x;
            /// <summary>
            /// Final direction leaving the medium
            /// </summary>
            public float3 w;
            /// <summary>
            /// A representative position where the path scatters. If X is the k-th position, then the pdf is Phi^k / (1 + Phi + ... +  Phi^N).
            /// </summary>
            public float3 X;
            /// <summary>
            /// The direction arriving to the represenative position
            /// </summary>
            public float3 W;
        }

        /// <summary>
        /// Settings of the medium ruling the experiment
        /// </summary>
        public class MediumSettings
        {
            /// <summary>
            /// Gets the density (extinction coefficient) of the medium
            /// </summary>
            public float Sigma;

            /// <summary>
            /// Gets the scattering albedo of the medium
            /// </summary>
            public float Phi;

            /// <summary>
            /// Gets the HG factor of the phase function
            /// </summary>
            public float G;
        }

        static float invertcdf(float g, float xi)
        {
            float t = (1.0f - g*g) / (1.0f - g + 2.0f * g * xi);
            return 0.5f * (1 + g * g - t * t) / g;
        }

        /// <summary>
        /// Implementation of the HG phase function.
        /// </summary>
        /// <param name="w">Incomming direction</param>
        /// <param name="g">Anisotropy factor</param>
        static float3 SamplePhase(float3 w, float g)
        {
            if (abs(g) < 0.001f)
                return randomDirection(-w);

            float phi = random() * 2 * pi;
            float cosTheta = invertcdf(g, random());
            float sinTheta = sqrt(max(0, 1.0f - cosTheta * cosTheta));

            createOrthoBasis(w, out float3 t0, out float3 t1);

            return sinTheta * sin(phi) * t0 + sinTheta * cos(phi) * t1 +
                cosTheta * w;
        }

        /// <summary>
        /// Distance needed to leave the sphere from x in direction d.
        /// </summary>
        static float DistanceToBoundary (float3 x, float3 d)
        {
            //float a = dot(d,d); <- 1 because d is normalized
            float b = 2 * dot(x, d);
            float c = dot(x, x) - 1;

            float Disc = b * b - 4 * c;

            if (Disc <= 0)
                return 0;

            // Assuming x is inside the sphere, only the positive root is needed (intersection forward w).
            return max(0, (-b + sqrt(Disc)) / 2); 
        }

        /// <summary>
        /// Performs a VPT in a sphere medium with specific settings starting at center and returns the summary of the path traced.
        /// </summary>
        static PathSummary GetVPTSampleInSphere (MediumSettings settings)
        {
            float3 x = float3(0, 0, 0);
            float3 w = float3(0, 0, 1);
            float3 X = x;
            float3 W = w;
            int N = 0;
            float accum = 0;
            float importance = 1;

            while (true)
            {
                importance *= settings.Phi;
                accum += importance;

                if (random() < importance / accum) // replace the representative by this one
                {
                    X = x;
                    W = w;
                }

                w = SamplePhase(w, settings.G);

                N++;

                float d = DistanceToBoundary(x, w);

                float t = settings.Sigma < 0.00001 ? 10000000 : -log(max(0.000000001f, 1 - random())) / settings.Sigma;

                if (t >= d || float.IsNaN(t) || float.IsInfinity(t))
                {
                    x += w * d;
                    return new PathSummary
                    {
                        N = N,
                        x = x,
                        w = w,
                        X = X,
                        W = W
                    };
                }
                x += w * t;
            }
        }

        static MediumSettings GenerateNewSettings()
        {
            float densityRnd = random();
            float scatterAlbedoRnd = random();
            float gRnd = random();

            float density = pow(densityRnd, 2) * 300; // densities varies from 0 to 300
            float scatterAlbedo = min(1, 1.000001f - pow(scatterAlbedoRnd, 6)); // transform the albedo testing set to vary really slow close to 1.
            float g = clamp(gRnd * 2 - 1, -0.999f, 0.999f); // avoid singular cases abs(g)=1

            return new MediumSettings
            {
                Sigma = density,
                Phi = scatterAlbedo,
                G = g
            };
        }

        static void Main(string[] args)
        {
            const int N = 1 << 22;

            StreamWriter writer = new StreamWriter("ScattersDataSet.ds");
            Console.WriteLine("Generating file...");

            for (int i=0; i<N; i++)
            {
                if (i % 1000 == 0)
                {
                    Console.Write("\r                                                ");
                    Console.Write("\rCompleted... " + (i * 100.0f / N).ToString("F2"));
                }

                var settings = GenerateNewSettings();
                var r = GetVPTSampleInSphere(settings);

                // code path variables in a compact way
                float3 zAxis = float3(0, 0, 1);
                float3 xAxis = abs(r.x.z) > 0.999 ? float3(1, 0, 0) : normalize(cross(r.x, float3(0, 0, 1)));
                float3 yAxis = cross(zAxis, xAxis);

                float3x3 normR = transpose(float3x3(xAxis, yAxis, zAxis));
                float3 normx = mul(r.x, normR);
                float3 normw = mul(r.w, normR);
                float3 normX = mul(r.X, normR);
                float3 normW = mul(r.W, normR);
                float3 B = float3(1, 0, 0);
                float3 T = cross(normx, B);
                float costheta = normx.z;
                float beta = dot(normw, T);
                float alpha = dot(normw, B);

                writer.WriteLine("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}",
                    settings.Sigma,
                    settings.G,
                    settings.Phi,
                    r.N,
                    costheta,
                    beta,
                    alpha,
                    normX.x,
                    normX.y,
                    normX.z,
                    normW.x,
                    normW.y,
                    normW.z
                    );
            }
            Console.WriteLine();

            writer.Close();
            Console.WriteLine("Done.");
        }
    }
}
