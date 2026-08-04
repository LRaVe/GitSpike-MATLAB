#include "mex.h"

#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

struct ProfilePoint
{
    double t;
    double s;
};

///////////////////////////////////////////////////////////////////////////
// Convert MATLAB vector -> std::vector<double>
///////////////////////////////////////////////////////////////////////////

static vector<double> getVector(const mxArray* A)
{
    size_t N = mxGetNumberOfElements(A);

    double* p = mxGetPr(A);

    return vector<double>(p,p+N);
}

///////////////////////////////////////////////////////////////////////////
// Equivalent de find(...,'last')
///////////////////////////////////////////////////////////////////////////

static int findLastLE(const vector<double>& v,double x)
{
    int idx=0;

    for(size_t i=0;i<v.size();i++)
    {
        if(v[i]<=x)
            idx=(int)i;
        else
            break;
    }

    return idx;
}

///////////////////////////////////////////////////////////////////////////
// auxiliary_delta (traduction fidèle Matlab)
///////////////////////////////////////////////////////////////////////////

static double auxiliary_delta(
        double spike,
        const vector<double>& own_train,
        const vector<double>& other_train,
        int idx,
        bool aux_idx)
{
    double delta=fabs(spike-other_train[0]);

    for(size_t i=1;i<other_train.size();i++)
        delta=min(delta,fabs(spike-other_train[i]));

    if(idx==0 && aux_idx)
    {
        delta=fabs(own_train[1]-other_train[0]);

        for(size_t i=1;i<other_train.size();i++)
            delta=min(delta,fabs(own_train[1]-other_train[i]));
    }

    if(idx==(int)own_train.size()-1 && aux_idx)
    {
        delta=fabs(own_train[own_train.size()-2]-other_train[0]);

        for(size_t i=1;i<other_train.size();i++)
            delta=min(delta,fabs(own_train[own_train.size()-2]-other_train[i]));
    }

    return delta;
}

///////////////////////////////////////////////////////////////////////////
// mexFunction
///////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs,
                 mxArray* plhs[],
                 int nrhs,
                 const mxArray* prhs[])
{
    if(nrhs!=10)
        mexErrMsgTxt("Ten inputs expected.");

    vector<double> spikes1=getVector(prhs[0]);
    vector<double> spikes2=getVector(prhs[1]);

    double t_min=mxGetScalar(prhs[2]);
    double t_max=mxGetScalar(prhs[3]);

    bool aux1_begin=(mxGetScalar(prhs[4])!=0);
    bool aux1_end=(mxGetScalar(prhs[5])!=0);
    bool aux2_begin=(mxGetScalar(prhs[6])!=0);
    bool aux2_end=(mxGetScalar(prhs[7])!=0);

    vector<double> Distances=getVector(prhs[8]);

    double threshold=mxGetScalar(prhs[9]);

    //////////////////////////////////////////////////////////////////////
    // SPIKE_distance_profile = cell(1,4)
    //////////////////////////////////////////////////////////////////////

    vector<ProfilePoint> SPIKE_distance_profile[4];

    //////////////////////////////////////////////////////////////////////
    // FIRST MAIN LOOP
    //////////////////////////////////////////////////////////////////////

    for(int idx1=0; idx1<(int)spikes1.size(); idx1++)
    {
        int idx2;

        if(spikes2.front()>spikes1[idx1])
            idx2=0;
        else if(spikes2.back()<=spikes1[idx1])
            idx2=(int)spikes2.size()-2;
        else
            idx2=findLastLE(spikes2,spikes1[idx1]);

        /////////////////////////////////////////////////////////////////
        // Train 2 contribution
        /////////////////////////////////////////////////////////////////

        double ISI_dist_2=
            spikes2[idx2+1]-spikes2[idx2];

        double delta_tp_2=
            auxiliary_delta(
                spikes2[idx2],
                spikes2,
                spikes1,
                idx2,
                aux2_begin);

        double delta_tf_2=
            auxiliary_delta(
                spikes2[idx2+1],
                spikes2,
                spikes1,
                idx2+1,
                aux2_end);

        double xp_2=
            spikes1[idx1]-spikes2[idx2];

        double xf_2=
            spikes2[idx2+1]-spikes1[idx1];

        double S_2=
            ((delta_tp_2*xf_2)
            +(delta_tf_2*xp_2))
            /ISI_dist_2;

        /////////////////////////////////////////////////////////////////
        // MATLAB : if idx_1>1
        /////////////////////////////////////////////////////////////////

        if(idx1>0)
        {
            double ISI_dist_1=
                spikes1[idx1]-spikes1[idx1-1];

            double S_1=
                auxiliary_delta(
                    spikes1[idx1],
                    spikes1,
                    spikes2,
                    idx1,
                    aux1_end);

            double meanISI=
                0.5*(ISI_dist_1+ISI_dist_2);

            if(Distances[0])
            {
                double S=
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*meanISI);

                SPIKE_distance_profile[0].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[1])
            {
                double S=
                    (S_1+S_2)
                    /(2.0*meanISI);

                SPIKE_distance_profile[1].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[2])
            {
                double denom=
                    max(meanISI,threshold);

                double S=
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*denom);

                SPIKE_distance_profile[2].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[3])
            {
                double denom=
                    max(meanISI,threshold);

                double S=
                    (S_1+S_2)
                    /(2.0*denom);

                SPIKE_distance_profile[3].push_back(
                    {spikes1[idx1],S});
            }
        }

        /////////////////////////////////////////////////////////////////
        // MATLAB : if idx_1 < length(spikes1)
        /////////////////////////////////////////////////////////////////

        if(idx1 < (int)spikes1.size()-1)
        {
            double ISI_dist_1 =
                spikes1[idx1+1]-spikes1[idx1];

            double S_1 =
                auxiliary_delta(
                    spikes1[idx1],
                    spikes1,
                    spikes2,
                    idx1,
                    aux1_begin);

            double meanISI =
                0.5*(ISI_dist_1+ISI_dist_2);

            if(Distances[0])
            {
                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*meanISI);

                SPIKE_distance_profile[0].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[1])
            {
                double S =
                    (S_1+S_2)
                    /(2.0*meanISI);

                SPIKE_distance_profile[1].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[2])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*denom);

                SPIKE_distance_profile[2].push_back(
                    {spikes1[idx1],S});
            }

            if(Distances[3])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    (S_1+S_2)
                    /(2.0*denom);

                SPIKE_distance_profile[3].push_back(
                    {spikes1[idx1],S});
            }
        }

    } // FIN de la boucle idx1

    //////////////////////////////////////////////////////////////////////
    // SECOND MAIN LOOP (idx2)
    //////////////////////////////////////////////////////////////////////

    for(int idx2=0; idx2<(int)spikes2.size(); idx2++)
    {
        int idx1;

        if(spikes1.front()>spikes2[idx2])
            idx1=0;
        else if(spikes1.back()<=spikes2[idx2])
            idx1=(int)spikes1.size()-2;
        else
            idx1=findLastLE(spikes1,spikes2[idx2]);

        /////////////////////////////////////////////////////////////////
        // Train 1 contribution
        /////////////////////////////////////////////////////////////////

        double ISI_dist_1 =
            spikes1[idx1+1]-spikes1[idx1];

        double delta_tp_1 =
            auxiliary_delta(
                spikes1[idx1],
                spikes1,
                spikes2,
                idx1,
                aux1_begin);

        double delta_tf_1 =
            auxiliary_delta(
                spikes1[idx1+1],
                spikes1,
                spikes2,
                idx1+1,
                aux1_end);

        double xp_1 =
            spikes2[idx2]-spikes1[idx1];

        double xf_1 =
            spikes1[idx1+1]-spikes2[idx2];

        double S_1 =
            ((delta_tp_1*xf_1)
            +(delta_tf_1*xp_1))
            /ISI_dist_1;

        /////////////////////////////////////////////////////////////////
        // MATLAB : if idx_2 > 1
        /////////////////////////////////////////////////////////////////

        if(idx2>0)
        {
            double ISI_dist_2 =
                spikes2[idx2]-spikes2[idx2-1];

            double S_2 =
                auxiliary_delta(
                    spikes2[idx2],
                    spikes2,
                    spikes1,
                    idx2,
                    aux2_end);

            double meanISI =
                0.5*(ISI_dist_1+ISI_dist_2);

            if(Distances[0])
            {
                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*meanISI);

                SPIKE_distance_profile[0].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[1])
            {
                double S =
                    (S_1+S_2)
                    /(2.0*meanISI);

                SPIKE_distance_profile[1].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[2])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*denom);

                SPIKE_distance_profile[2].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[3])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    (S_1+S_2)
                    /(2.0*denom);

                SPIKE_distance_profile[3].push_back(
                    {spikes2[idx2],S});
            }
        }

        /////////////////////////////////////////////////////////////////
        // MATLAB : if idx_2 < length(spikes2)
        /////////////////////////////////////////////////////////////////

        if(idx2<(int)spikes2.size()-1)
        {
            double ISI_dist_2 =
                spikes2[idx2+1]-spikes2[idx2];

            double S_2 =
                auxiliary_delta(
                    spikes2[idx2],
                    spikes2,
                    spikes1,
                    idx2,
                    aux2_begin);

            double meanISI =
                0.5*(ISI_dist_1+ISI_dist_2);

            if(Distances[0])
            {
                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*meanISI);

                SPIKE_distance_profile[0].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[1])
            {
                double S =
                    (S_1+S_2)
                    /(2.0*meanISI);

                SPIKE_distance_profile[1].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[2])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    ((S_1*ISI_dist_2)
                    +(S_2*ISI_dist_1))
                    /(2.0*meanISI*denom);

                SPIKE_distance_profile[2].push_back(
                    {spikes2[idx2],S});
            }

            if(Distances[3])
            {
                double denom =
                    max(meanISI,threshold);

                double S =
                    (S_1+S_2)
                    /(2.0*denom);

                SPIKE_distance_profile[3].push_back(
                    {spikes2[idx2],S});
            }
        }

    } // fin de la boucle idx2

    //////////////////////////////////////////////////////////////////////
    // CONVERT PROFILE (traduction directe Matlab)
    //////////////////////////////////////////////////////////////////////

    vector< vector<ProfilePoint> > profile_mat(4);

    for(int k=0;k<4;k++)
    {
        if(!Distances[k])
            continue;

        vector<ProfilePoint> profile =
            SPIKE_distance_profile[k];

        /////////////////////////////////////////////////////////////////
        // sortrows(profile,1)
        /////////////////////////////////////////////////////////////////

        stable_sort(profile.begin(),
             profile.end(),
             [](const ProfilePoint& a,
                const ProfilePoint& b)
             {
                 return a.t < b.t;
             });

        /////////////////////////////////////////////////////////////////
        // keep only points inside interval
        /////////////////////////////////////////////////////////////////

        for(size_t i=0;i<profile.size();i++)
        {
            /////////////////////////////////////////////////////////////
            // point before t_min
            /////////////////////////////////////////////////////////////

            if(profile[i].t < t_min)
            {
                size_t idx=0;

                while(idx<profile.size() &&
                      profile[idx].t<t_min)
                    idx++;

                if(idx<profile.size())
                {
                    profile[i].s =
                        profile[i].s +
                        ((profile[idx].s-profile[i].s)
                        /(profile[idx].t-profile[i].t))
                        *(t_min-profile[i].t);
                }

                profile[i].t=t_min;
            }

            /////////////////////////////////////////////////////////////
            // point after t_max
            /////////////////////////////////////////////////////////////

            else if(profile[i].t>t_max)
            {
                int idx=(int)profile.size()-1;

                while(idx>=0 &&
                      profile[idx].t>t_max)
                    idx--;

                if(idx>=0)
                {
                    profile[i].s =
                        profile[idx].s +
                        ((profile[i].s-profile[idx].s)
                        /(profile[i].t-profile[idx].t))
                        *(t_max-profile[idx].t);
                }

                profile[i].t=t_max;
            }
        }

        /////////////////////////////////////////////////////////////////
        // unique(...,'rows','stable')
        /////////////////////////////////////////////////////////////////

        vector<ProfilePoint> uniqueProfile;

        for(size_t i=0;i<profile.size();i++)
        {
            bool alreadySeen=false;
        
            for(size_t j=0;j<uniqueProfile.size();j++)
            {
                if(profile[i].t==uniqueProfile[j].t &&
                   profile[i].s==uniqueProfile[j].s)
                {
                    alreadySeen=true;
                    break;
                }
            }
        
            if(!alreadySeen)
                uniqueProfile.push_back(profile[i]);
        }

        /////////////////////////////////////////////////////////////////
        // sortrows(profile,1)
        /////////////////////////////////////////////////////////////////

        stable_sort(uniqueProfile.begin(),
             uniqueProfile.end(),
             [](const ProfilePoint& a,
                const ProfilePoint& b)
             {
                 return a.t < b.t;
             });

        profile_mat[k]=uniqueProfile;
    }


    //////////////////////////////////////////////////////////////////////
    // FINAL DISTANCE (trapz)
    //////////////////////////////////////////////////////////////////////

    double SPIKE_distance_2x2[4]={0,0,0,0};

    for(int k=0;k<4;k++)
    {
        if(!Distances[k])
            continue;

        const auto& profile=profile_mat[k];

        double area=0.0;

        for(size_t i=0;i+1<profile.size();i++)
        {
            area +=
                (profile[i+1].t-profile[i].t)*
                (profile[i].s+profile[i+1].s)/2.0;
        }

        SPIKE_distance_2x2[k]=
            area/(t_max-t_min);
    }

    //////////////////////////////////////////////////////////////////////
    // OUTPUT 1 : SPIKE_distance_2x2
    //////////////////////////////////////////////////////////////////////

    plhs[0]=mxCreateDoubleMatrix(1,4,mxREAL);

    double* out=
        mxGetPr(plhs[0]);

    for(int k=0;k<4;k++)
        out[k]=SPIKE_distance_2x2[k];

    //////////////////////////////////////////////////////////////////////
    // OUTPUT 2 : profile_mat
    //////////////////////////////////////////////////////////////////////

    plhs[1]=mxCreateCellMatrix(1,4);

    for(int k=0;k<4;k++)
    {
        if(!Distances[k] || profile_mat[k].empty())
            continue;

        mwSize rows=profile_mat[k].size();

        mxArray* mat=
            mxCreateDoubleMatrix(rows,2,mxREAL);

        double* ptr=
            mxGetPr(mat);

        for(mwSize i=0;i<rows;i++)
        {
            ptr[i]          = profile_mat[k][i].t;
            ptr[i+rows]     = profile_mat[k][i].s;
        }

        mxSetCell(plhs[1],k,mat);
    }

} // ===== FIN mexFunction =====