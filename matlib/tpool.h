#ifndef THREAD_POOL_H_
#define THREAD_POOL_H_

#include <array>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>

template <unsigned> class thread_pool;

namespace internal {

template <unsigned I> void thread_pool_worker(thread_pool<I> *tp) {
  while (tp->m_running) {
    std::function<void()> task;
    {
      std::unique_lock<std::mutex> lock(tp->m_mutex);
      tp->m_condition.wait(
          lock, [tp] { return !tp->m_running || !tp->m_jobs.empty(); });
      if(!tp->m_running && tp->m_jobs.empty())
        return;
      task = std::move(tp->m_jobs.front());
      tp->m_jobs.pop();
    }
    task();
  }
}

} // namespace internal

template <unsigned I> class thread_pool {
public:
  inline thread_pool() : m_running(true) {
    for (int i = 0; i < I; i++) {
      m_threads[i] = std::thread(internal::thread_pool_worker<I>, this);
    }
  }
  inline ~thread_pool() {
    m_running = false;
    m_condition.notify_all();
    for (int i = 0; i < I; i++)
      m_threads[i].join();
  }
  template <typename F, typename... Args>
  auto operator()(F &&f, Args &&... args)
      -> std::future<typename std::result_of<F(Args...)>::type> {
    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    std::future<return_type> res = task->get_future();
    {
      std::unique_lock<std::mutex> lock(m_mutex);

      m_jobs.push([task]() { (*task)(); });
    }
    m_condition.notify_one();
    return res;
  }

private:
  template <unsigned S>
  friend void internal::thread_pool_worker(thread_pool<S> *tp);
  std::array<std::thread, I> m_threads;
  std::queue<std::function<void()>> m_jobs;
  std::mutex m_mutex;
  bool m_running;
  std::condition_variable m_condition;
};

#endif
