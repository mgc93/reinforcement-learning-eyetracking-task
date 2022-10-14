function [isi_t1a_1, isi_t2a_1,isi_t1a_2, isi_t2a_2, eff_save_1a_1, eff_save_2a_1,eff_save_alla_1, eff_save_1a_2, eff_save_2a_2,eff_save_alla_2] = createTrialTiming(nsim, lambda, T, min_iti)

% to do list

% add an impulse function for the response (distributed truncated exgaussian? or
% uniformly distributed with min = 2 and max = 3)



%looking to see what jitter does to efficiency in TR's of 2.6s

%addpath /Users/miruna/Documents/MATLAB/Project learning/
%edit the above path to point towards the directory with spm_hrf in it.


% goal is to have a design with a stimulus followed by feedback.


% How many stimuli can be fit into one 1200s long run? (about 20 minutes)



%% design with uniform jitter

% %Uniform distribution between 0 & 1: randomly select numbers between 0 and 1.  Mean of
% %U(0,1) is 0.5, mean of U(1,4) is (4+1)/2, mean of U(2,6) is (6+2)/2=8
% 
% 
% 
% % keep the first fixation 4s on average and the second
% %second 4s on average. If the minimum fixation is 2s long, the uniform
% %distribution range for an average duration of 4s?
% 
% 
% %U(2,6) = U(0,1)*4 + 2
% 
% 
% %the duration of the stimulus is 3s long, so add
% %that to the Uniform distribution.
% 
% %create each trial in high resolution, convolve and then downsample to
% %a TR of TR second
% %set up the hrf info
% TR=2.6; % as in Konovalov and Krajbich (2018)
% 
% t=0:0.2:d; % 0.2 the resolution should not be = TR, but much higher
% % SPM TR/16, here use 0.2 = TR/13
% hrf_25=spm_hrf(0.2);
% 
% % create the time for stimulus(t1) and feedback(t2)
% t1b(1) = 1;
% for i = 1:105
% 
%     t2b(i)   = t1b(i)+3+rand(1)*4+2;
%     % when the last stimulus exceeds the fixed duration of the experiment
%     % do not create additional stimuli
%     if(t2b(i)>d)
%         break;
%     end
%     t1b(i+1) = t2b(i)+2+rand(1)*4+2;
% 
% 
% end
% 
% r1b = zeros(1, d*5+1);  %I'm assuming time resoultion of .2 s, so this corresponds to 875 2.6s TRs
% for i = 1:length(t1b)
%     r1b(t1b(i)<=t & t<=(t1b(i)+3))=1;  %add 3 for a 3s trial
% end
% 
% r1b = conv(hrf_25, r1b);
% r1b = r1b(1:5*TR:(d*5)); % sample it down to the TR resolution
% 
% t_tr = t(1:5*TR:(d*5));  %this is time in seconds (for plotting purposes)
% 
% 
% 
% r2b = zeros(1, (d*5+1));
% 
% for i = 1:length(t2b)
%   r2b(t2b(i)<=t & t<=(t2b(i)+3))=1;  %add 3 for a 3s trial
% end
% 
% r2b = conv(hrf_25, r2b);
% r2b = r2b(1:5*TR:(d*5));
% 
% % mean centering, or put a column of 1s in the design matrix - intercept
% r1b = r1b -mean(r1b);
% r2b = r2b-mean(r2b);
% 
% 
% 
% figure();
% 
% plot(t_tr, r1b, 'linewidth', 2)
% hold on
% plot(t_tr, r2b, 'g', 'linewidth', 2)
% hold off
% xlabel('Time (s)')
% title('Jitter uniform between 2 and 6 s', 'Fontsize', 14)
% 
% %compare the efficiences of the 2 designs
% Xb=[r1b', r2b'];
% 
% c1=[1 0]; % contrast for stimulus alone
% c2=[0 1]; % contrast for response alone
% eff1b=1./(c1*inv(Xb'*Xb)*c1');
% 
% eff2b=1./(c2*inv(Xb'*Xb)*c2');
% 
% eff_allb=2./(c1*inv(Xb'*Xb)*c1'+c2*inv(Xb'*Xb)*c2');
% 
% % fprintf('c1:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff1b, eff1);
% % fprintf('c2:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff2b, eff2);
% %
% % fprintf('all:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff_allb, eff_all);
% corr(r1b', r2b');
% corr(r1b', r2b');
% 
% 
% %nsim=100;
% eff_save_1b=zeros(nsim, 1);
% eff_save_2b=zeros(nsim,1);
% 
% %  save the ISI information for each simulation to be able
% % to access it later
% isi_sim_t1b = [];
% isi_sim_t2b = [];
% 
% for j=1:nsim
%     
%     % create the time for stimulus(t1) and feedback(t2)
%     t1b(1) = 1;
%     for i = 1:105
%         
%         t2b(i)   = t1b(i)+3+rand(1)*4+2;
%         
%     % when the last stimulus exceeds the fixed duration of the experiment
%     % do not create additional stimuli
%     if(t2b(i)>d)
%         break;
%     end
%         t1b(i+1) = t2b(i)+2+rand(1)*4+2;
%         
%     end
%     
%     r1b = zeros(1, d*5+1);  %I'm assuming time resoultion of .2 s, so this corresponds to 875 2.6s TRs
%     for i = 1:length(t1b)
%         r1b(t1b(i)<=t & t<=(t1b(i)+3))=1;  %add 3 for a 3s trial
%     end
%     
%     r1b = conv(hrf_25, r1b);
%     r1b = r1b(1:5*TR:(d*5)); % sample it down to the TR resolution
%     
%     t_tr = t(1:5*TR:(d*5));  %this is time in seconds (for plotting purposes)
%     
%     
%     
%     r2b = zeros(1, (d*5+1));
%     
%     for i = 1:length(t2b)
%         r2b(t2b(i)<=t & t<=(t2b(i)+3))=1;  %add 3 for a 3s trial
%     end
%     
%     r2b = conv(hrf_25, r2b);
%     r2b = r2b(1:5*TR:(d*5));
%     
%     % mean centering, or put a column of 1s in the design matrix - intercept
%     r1b = r1b-mean(r1b);
%     r2b = r2b-mean(r2b);
%     
%     
%     %compare the efficiences of the 2 designs
%     Xb=[r1b', r2b'];
%     c1=[1 0]; % contrast for stimulus alone
%     c2=[0 1]; % contrast for response alone
%     eff_save_1b(j,1)=1./(c1*inv(Xb'*Xb)*c1');
%     eff_save_2b(j,1)=1./(c2*inv(Xb'*Xb)*c2');
%     eff_save_allb(j,1) = 2./(c1*inv(Xb'*Xb)*c1'+c2*inv(Xb'*Xb)*c2');
%     
%     isi_sim_t1b = [isi_sim_t1b, t1b'];
%     isi_sim_t2b = [isi_sim_t2b, t2b'];
%     
% end
% 
% % subplot(1,1,1)
% % plot(eff_save_1b, eff_save_2b, '.')

%% design with truncated exponential jitter

% % %% FYI, here's how you properly code up a truncated exponential
% %
% n = 1000000;
% lambda = 6.12;    % This is the parameter on your exponential
% T = 8;            % This is the largest duration allowed
% 
% % Here's the mean of the truncated exponential
% mean_te = lambda - T*(exp(T/lambda)-1)^(-1);
% 
% R_1 = rand(n,1)*(1-exp(-T/lambda));
% R_2 = rand(n,1)*(1-exp(-T/lambda));
% rand_isi_1 = -log(1-R_1)*lambda;
% rand_isi_2 = -log(1-R_2)*lambda;
% 
% %mean(rand_isi)
% 
% %create each trial in high resolution, convolve and then downsample to
% %a TR of TR second
% %set up the hrf info
% TR=2.6; % as in Konovalov and Krajbich (2018)
% 
% t=0:0.2:d; % 0.2 the resolution should not be = TR, but much higher
% % SPM TR/16, here use 0.2 = TR/13
% hrf_25=spm_hrf(0.2);
% 
% % create the time for stimulus(t1) and feedback(t2)
% t1a(1) = 1;
% for i = 1:105
%     
%     t2a(i)=t1a(i)+3+rand_isi_2(i);
%     
%     % when the last stimulus exceeds the fixed duration of the experiment
%     % do not create additional stimuli
%     if(t2a(i)>d)
%         break;
%     end
%     t1a(i+1) = t2a(i)+2+rand_isi_1(i);
%     
% end
% 
% r1a = zeros(1, (d*5+1));  %I'm assuming time resoultion of .2 s, so this corresponds to 875 2.6s TRs
% for i=1:length(t1a)
%     r1a(t1a(i)<=t & t<=(t1a(i)+3))=1;  %add 3 for a 3s trial
% end
% 
% r1a = conv(hrf_25, r1a);
% r1a = r1a(1:5*TR:(d*5)); % sample it down to the TR resolution
% 
% t_tr=t(1:5*TR:(d*5));  %this is time in seconds (for plotting purposes)
% 
% 
% 
% r2a=zeros(1, (d*5+1));
% 
% for i=1:length(t2a)
%     r2a(t2a(i)<=t & t<=(t2a(i)+3))=1;  %add 3 for a 3s trial
% end
% 
% r2a = conv(hrf_25, r2a);
% r2a = r2a(1:5*TR:(d*5+1));
% 
% % mean centering, or put a column of 1s in the design matrix - intercept
% r1a = r1a-mean(r1a);
% r2a = r2a-mean(r2a);
% 
% 
% 
% figure();
% 
% plot(t_tr, r1a, 'linewidth', 2)
% hold on
% plot(t_tr, r2a, 'g', 'linewidth', 2)
% hold off
% xlabel('Time (s)')
% title('Jitter truncated exponential with T = 8, lambda = 6.12', 'Fontsize', 14)
% 
% %compare the efficiences of the 2 designs
% Xa=[r1a', r2a'];
% 
% c1=[1 0]; % contrast for stimulus alone
% c2=[0 1]; % contrast for response alone
% eff1a=1./(c1*inv(Xa'*Xa)*c1');
% 
% eff2a=1./(c2*inv(Xa'*Xa)*c2');
% 
% eff_alla=2./(c1*inv(Xa'*Xa)*c1'+c2*inv(Xa'*Xa)*c2');
% 
% % fprintf('c1:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff1b, eff1);
% % fprintf('c2:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff2b, eff2);
% %
% % fprintf('all:  Jittered efficiency=%g     Fixed efficiency=%g \n', eff_allb, eff_all);
% corr(r1a', r2a');
% corr(r1a', r2a');



nsim = 1000;
% d = 1200;
d = 1500;
min_iti = 2;


%nsim=100;
eff_save_1a=zeros(nsim,1);
eff_save_2a=zeros(nsim,1);

%  save the ISI information for each simulation to be able
% to access it later
isi_sim_t1a = [];
isi_sim_t2a = [];
for j=1:nsim
    % %% FYI, here's how you properly code up a truncated exponential
%
n = 1000000;
lambda = 2;    % This is the parameter on your exponential
T = 4;            % This is the largest duration allowed

% Here's the mean of the truncated exponential
mean_te = lambda - T*(exp(T/lambda)-1)^(-1);

R_1 = rand(n,1)*(1-exp(-T/lambda));
R_2 = rand(n,1)*(1-exp(-T/lambda));
rand_isi_1 = -log(1-R_1)*lambda;
rand_isi_2 = -log(1-R_2)*lambda;

mean_isi_1 = mean(rand_isi_1);
mean_isi_2 = mean(rand_isi_2);

%create each trial in high resolution, convolve and then downsample to
%a TR of TR second
%set up the hrf info
TR=2.6; % as in Konovalov and Krajbich (2018)

t=0:0.2:d; % 0.2 the resolution should not be = TR, but much higher
% SPM TR/16, here use 0.2 = TR/13
hrf_25=spm_hrf(0.2);

% create the time for stimulus(t1) and feedback(t2)
t1a(1) = 1;
for i = 1:105
    
    t2a(i)=t1a(i)+min_iti+3+rand_isi_2(i); % added 3s for the choice screen
    % added another 2 for the minimum
    
    % when the last stimulus exceeds the fixed duration of the experiment
    % do not create additional stimuli
    if(t2a(i)>d)
        break;
    end
    t1a(i+1) = t2a(i)+min_iti+2+rand_isi_1(i); % added 2s for the feedback screen
    % added another 2 for the minimum
end

    t1a = t1a(1,[1:105]);


r1a = zeros(1, (d*5+1));  %I'm assuming time resoultion of .2 s, so this corresponds to 875 2.6s TRs
for i=1:length(t1a)
    r1a(t1a(i)<=t & t<=(t1a(i)+3))=1;  %add 3 for a 3s trial
end

r1a = conv(hrf_25, r1a);
r1a = r1a(1:5*TR:(d*5)); % sample it down to the TR resolution

t_tr=t(1:5*TR:(d*5));  %this is time in seconds (for plotting purposes)



r2a=zeros(1, (d*5+1));

for i=1:length(t2a)
    r2a(t2a(i)<=t & t<=(t2a(i)+3))=1;  %add 3 for a 3s trial
end

r2a = conv(hrf_25, r2a);
r2a = r2a(1:5*TR:(d*5));

% mean centering, or put a column of 1s in the design matrix - intercept
r1a = r1a-mean(r1a);
r2a = r2a-mean(r2a);
    
    Xa = [r1a', r2a'];
    c1=[1 0]; % contrast for stimulus alone
    c2=[0 1]; % contrast for response alone
    eff_save_1a(j,1)=1./(c1*inv(Xa'*Xa)*c1');
    eff_save_2a(j,1)=1./(c2*inv(Xa'*Xa)*c2');
    eff_save_alla(j,1) = 2./(c1*inv(Xa'*Xa)*c1'+c2*inv(Xa'*Xa)*c2');
    
    isi_sim_t1a = [isi_sim_t1a, t1a'];
    isi_sim_t2a = [isi_sim_t2a, t2a'];
    
end

% subplot(1,1,1)
% plot(eff_save_1, eff_save_2, '.')

% cut the length of isi_t1a to be the same as isi_t2a


% select 3 simulation with maximum efficiency

%run 1
[m1, ind1] = max(eff_save_alla);
isi_t1a_1 = isi_sim_t1a(:,ind1);
isi_t2a_1 = isi_sim_t2a(:,ind1);
eff_save_1a_1 = eff_save_1a(ind1,1);
eff_save_2a_1 = eff_save_2a(ind1,1);
eff_save_alla_1 = eff_save_alla(ind1,1);

isi_sim_t1a(:,ind1) = [];
isi_sim_t2a(:,ind1) = [];
eff_save_1a(ind1,:) = [];
eff_save_2a(ind1,:) = [];
eff_save_alla(ind1,:) = [];

[m2, ind2] = max(eff_save_alla);
isi_t1a_2 = isi_sim_t1a(:,ind2);
isi_t2a_2 = isi_sim_t2a(:,ind2);
eff_save_1a_2 = eff_save_1a(ind2,1);
eff_save_2a_2 = eff_save_2a(ind2,1);
eff_save_alla_2 = eff_save_alla(ind2,1);


% figure();
% plot(eff_save_1a, eff_save_2a, 'b.');
% xlabel('Efficiency for contrast for stimulus');
% ylabel('Efficiency for contrast for feedback');

%% compare both types

% figure();
% plot(eff_save_1a, eff_save_2a, 'b.');
% xlabel('Efficiency for contrast for stimulus');
% ylabel('Efficiency for contrast for feedback');
% hold on
% plot(eff_save_1b, eff_save_2b, 'r.');
% xlabel('Efficiency for contrast for stimulus');
% ylabel('Efficiency for contrast for feedback');
% hold off





%% compare different parameters for the truncated exponential


% figure();
% plot(eff_1a, eff_2a, 'b.');
% xlabel('Efficiency for contrast for stimulus');
% ylabel('Efficiency for contrast for feedback');
% 
% hold on
% plot(eff_1c, eff_2c, 'r.');
% xlabel('Efficiency for contrast for stimulus');
% ylabel('Efficiency for contrast for feedback');
% title('Truncated exponential timing');
% hold off
% legend('lambda = 4, T = 8','lambda = 4, T = 6')


%% plot the times for each stimulus
figure()
subplot(4,1,1)
plot(1:104,isi_t1a_1([2:105],1)-isi_t1a_1([1:104],1),'b')
hold on
plot(1:104,isi_t2a_1([2:105],1)-isi_t2a_1([1:104],1),'r')
hold off
title(['Trials length for stimulus 1 and 2 (set 1)']);

subplot(4,1,2)
plot(1:104,isi_t1a_2([2:105],1)-isi_t1a_2([1:104],1),'b')
hold on
plot(1:104,isi_t2a_2([2:105],1)-isi_t2a_2([1:104],1),'r')
hold off
title(['Trials length for stimulus 1 and 2 (set 2)']);

subplot(4,1,3)
plot(1:105,isi_t2a_1(:,1)-isi_t1a_1(:,1),'k')
title(['Time between stimuli (set 1)']);
subplot(4,1,4)
plot(1:105,isi_t2a_2(:,1)-isi_t1a_2(:,1),'k')
title(['Time between stimuli (set 2)']);

end




