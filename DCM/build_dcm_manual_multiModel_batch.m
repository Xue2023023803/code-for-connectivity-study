%% build_dcm_manual_multiModel_batch.m
% 为每个被试提交一个 batch 作业，
% 每个作业在独立 MATLAB 进程中运行 build_dcm_manual_multiModel_spm25({subj})

% 被试列表（和主函数里默认的一致，也可以在这里单独改）
subjects = {'sub-01','sub-02','sub-08'};

% 使用本机 local cluster
c = parcluster('local');

% 根据 CPU 核数和被试数给一个推荐的最大并行作业数
maxWorkers = min(numel(subjects), feature('numcores'));
fprintf('推荐 local.NumWorkers 不超过 %d（<= 被试数且 <= CPU 核心数）。\n', maxWorkers);

% 如果 MATLAB 版本允许，也可以尝试在代码里设置（否则会报 warning，可以忽略）
try
    c.NumWorkers = maxWorkers;
catch
    warning('无法在代码中修改 c.NumWorkers，请在 Parallel Preferences 里把 local 的 NumWorkers 设为 <= %d。', maxWorkers);
end

% 提交 job（每个被试一个 job），每个 job 内部不再开 pool（Pool=0）
jobs = cell(numel(subjects), 1);
for i = 1:numel(subjects)
    subj = subjects{i};
    jobs{i} = batch(c, @build_dcm_manual_multiModel_spm25, 0, {{subj}}, 'Pool', 0);
    fprintf('已提交 job #%d -> 被试 %s\n', i, subj);
end

% 把 job 对象丢到 base workspace，方便你之后查状态
assignin('base', 'dcm_jobs', jobs);
fprintf('所有 batch job 已提交，可用 dcm_jobs{i}.State 查看状态。\n');
