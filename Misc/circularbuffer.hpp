/*
Header file for QuantumMechanics::CircularBuffer:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/

#ifndef _CIRCULARBUFFER_H_
#define _CIRCULARBUFFER_H_

#include <algorithm> // for std::min
#include <atomic>

template <size_t Capacity>
class CircularBuffer
{
private:
	std::atomic<size_t> write_begin, write_end, read_begin;
	std::atomic<size_t> unread_size;
	static const size_t data_capacity = Capacity;
	std::unique_ptr<char[]> data;
	std::atomic<bool> write_lock;
	std::atomic<bool> read_lock;

public:
	//Mutex behavior
	void lock_writes()
	{
		while (bool unlocked = false; !write_lock.compare_exchange_weak(unlocked, true)); // Spin-lock gaurding a write lock.
	}
	void unlock_writes()
	{
		write_lock = false;
	}
	void lock_reads()
	{
		while (bool unlocked = false; !read_lock.compare_exchange_weak(unlocked, true)); // Spin-lock gaurding a write lock.
	}
	void unlock_reads()
	{
		read_lock = false;
	}
	void wait_for_writes()
	{
		while (copy.write_begin.load() != copy.write_end.load()); // Spin-lock gaurding the current writes
	}
	void lock()
	{
		lock_writes();
		lock_reads();
	}
	void unlock()
	{
		unlock_writes();
		unlock_reads();
	}

	// Default constructor
	CircularBuffer() : 
		write_begin(0),
		write_end(0),
		read_begin(0),
		unread_size(0),
		data(new char[data_capacity])
	{ }

	// Default destructor
	virtual ~CircularBuffer();

	// Copy constructor
	CircularBuffer(const CircularBuffer& copy) :
		data(new char[data_capacity])
	{
		copy.lock();
		write_begin = copy.write_begin;
		write_end = copy.write_end;
		read_begin = copy.read_begin;
		unread_size = copy.unread_size;
		std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
		copy.unlock();
	}

	// Copy assignment operator
	CircularBuffer& operator=(const CircularBuffer& other)
	{
		if (this != &copy)
		{
			copy.lock();
			write_begin = copy.write_begin;
			write_end = copy.write_end;
			read_begin = copy.read_begin;
			unread_size = copy.unread_size;
			std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
			copy.unlock();
		}

		return *this;
	}
	
	// Move constructor
	CircularBuffer(CircularBuffer&& temp)
	{
		copy.lock();
		write_begin = std::move(copy.write_begin);
		write_end = std::move(copy.write_end);
		read_begin = std::move(copy.read_begin);
		unread_size = std::move(copy.unread_size);
		data = std::move(copy.data);
		copy.unlock();
	}

	// Move assignment operator
	CircularBuffer& operator=(CircularBuffer&& other)
	{
		if (this != &copy)
		{
			copy.lock();
			write_begin = std::move(copy.write_begin);
			write_end = std::move(copy.write_end);
			read_begin = std::move(copy.read_begin);
			unread_size = std::move(copy.unread_size);
			data = std::move(copy.data);
			copy.unlock();
		}

		return *this;
	}

	size_t size() const { return unread_size.load(); }
	size_t capacity() const { return data_capacity; }

protected:
	// Return number of bytes written.
	size_t write(const char *write_data, size_t bytes)
	{
		if (bytes == 0) return 0;

		lock_writes();

		size_t bytes_to_write = std::min(bytes, data_capacity - unread_size);

		// Write in a single step
		if (bytes_to_write <= data_capacity - write_end)
		{
			size_t current_write_end = write_end.fetch_add((write_end.load() + bytes_to_write == data_capacity) ? bytes_to_write : -write_end.load());
			unread_size += bytes_to_write;
			unlock_writes();
			memcpy(data + current_write_end, write_data, bytes_to_write);
		}
		// Write in two steps
		else
		{
			size_t current_write_end = write_end.load();
			size_t size_1 = data_capacity - current_write_end.load();
			size_t size_2 = bytes_to_write - size_1;
			write_end = size_2;
			unread_size += bytes_to_write;
			unlock_writes();
			memcpy(data + current_write_end, write_data, size_1);
			memcpy(data, write_data + size_1, size_2);
		}

		return bytes_to_write;
	}

	// Return number of bytes read.
	size_t read(char *read_data, size_t bytes)
	{
		if (bytes == 0) return 0;

		lock_reads();

		size_t bytes_to_read = std::min(bytes, unread_size);

		// Read in a single step
		if (bytes_to_read <= data_capacity - current_read_begin.load())
		{
			size_t current_read_begin = read_begin.fetch_add((read_begin.load() + bytes_to_write == data_capacity) ? bytes_to_write : -read_begin.load());
			unread_size -= bytes_to_write;
			unlock_reads();
			memcpy(read_data, data + current_read_begin, bytes_to_read);
		}
		// Read in two steps
		else
		{
			size_t current_read_begin = read_begin.load();
			size_t size_1 = data_capacity - current_read_begin;
			size_t size_2 = bytes_to_read - size_1;
			read_begin = size_2;
			unread_size -= bytes_to_write;
			unlock_reads();
			memcpy(read_data, data + current_read_begin, size_1);
			memcpy(read_data + size_1, data, size_2);
		}

		return bytes_to_read;
	}
};

#endif //namespace _CIRCULARBUFFER_H_
