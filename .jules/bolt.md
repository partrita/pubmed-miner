## 2025-05-22 - [Parallelize citation fetching]
**Learning:** Sequential IO-bound tasks like fetching citations from external APIs can be a significant bottleneck. Using `ThreadPoolExecutor` for parallelization dramatically improves throughput.
**Action:** Always look for batch operations that can be parallelized when dealing with external API calls.

## 2025-05-22 - [Thread-safe rate limiting]
**Learning:** Simple rate limiting logic using `time.sleep()` is not sufficient when multiple threads are involved. Threads can simultaneously check the same condition and bypass the limit.
**Action:** Use a lock to reserve time slots for each request to ensure rate limits are strictly enforced across all threads.
