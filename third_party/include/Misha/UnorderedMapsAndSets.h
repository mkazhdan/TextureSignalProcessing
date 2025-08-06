/*
Copyright (c) 2025, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef UNORDERED_MAPS_AND_SETS
#define UNORDERED_MAPS_AND_SETS

#include <unordered_map>
#include <unordered_set>

namespace MishaK
{
	// If the Key type does not have a Hasher subclass, depend on std::hash to work its magic
	template< typename Key , typename Value , typename = void >
	struct _UnorderedMap
	{
		using Type = std::unordered_map< Key , Value >;
	};

	template< typename Key , typename = void >
	struct _UnorderedSet
	{
		using Type = std::unordered_set< Key >;
	};

	// Otherwise, use the Key::Hasher subclass
	template< typename Key , typename Value >
	struct _UnorderedMap< Key , Value , std::void_t< typename Key::Hasher > >
	{
		using Type = std::unordered_map< Key , Value , typename Key::Hasher >;
	};

	template< typename Key >
	struct _UnorderedSet< Key , std::void_t< typename Key::Hasher > >
	{
		using Type = std::unordered_set< Key , typename Key::Hasher >;
	};


	template< typename Key , typename Value >
	using UnorderedMap = typename _UnorderedMap< Key , Value >::Type;

	template< typename Key >
	using UnorderedSet = typename _UnorderedSet< Key >::Type;

};

#endif // UNORDERED_MAPS_AND_SETS
