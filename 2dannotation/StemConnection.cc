#include "StemConnection.h"

#include <cassert>

namespace annotate
{
	StemConnection::StemConnection()
	{
		mpStructure = NULL;
		meConnection = Stem::eUNDEFINED_CONNECTION;
	}

	StemConnection::StemConnection(
		const Stem& aStem,
		const Stem::enConnection& aeConnect)
	{
		mpStructure = &aStem;
		meConnection = aeConnect;
	}

	mccore::ResId StemConnection::nextId() const
		throw(mccore::NoSuchElementException)
	{
		mccore::ResId id;
		const Stem* pStem = dynamic_cast<const Stem*>(mpStructure);
		assert(NULL != pStem);
		if(NULL != pStem && meConnection != Stem::eUNDEFINED_CONNECTION)
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				id = pStem->basePairs().front().rResId;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				id = pStem->basePairs().front().fResId;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				id = pStem->basePairs().back().rResId;
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				id = pStem->basePairs().back().fResId;
				break;
			default:
				std::string strMsg("Invalid connections have no nextId");
				throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
				break;
			}
		}
		else
		{
			std::string strMsg("Invalid connections have no nextId");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}

		return id;
	}

	BasePair StemConnection::getPair() const
		throw(mccore::NoSuchElementException)
	{
		const BasePair* pBasePair = NULL;
		const Stem* pStem = dynamic_cast<const Stem*>(mpStructure);
		assert(NULL != pStem);
		if(NULL != pStem && meConnection != Stem::eUNDEFINED_CONNECTION)
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				pBasePair = &pStem->basePairs().front();
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
			case Stem::eSECOND_STRAND_BACK_PAIR:
				pBasePair = &pStem->basePairs().back();
				break;
			default:
				pBasePair = NULL;
			}
		}

		if(NULL == pBasePair)
		{
			std::string strMsg("Invalid connections have no connecting pairs");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		return *pBasePair;
	}
/*
	bool StemConnection::connects(const StemConnection& aConnection) const
	{
		bool bConnects = false;
		if(isValid() && aConnection.isValid())
		{
			if(mpStructure == aConnection.mpStructure)
			{
				bConnects = (nextId() == aConnection.getResidue());
			}
		}
		return bConnects;
	}*/

	mccore::ResId StemConnection::getResidue() const
	{
		mccore::ResId resId;
		const Stem* pStem = dynamic_cast<const Stem*>(mpStructure);
		assert(NULL != pStem);
		if(NULL != pStem)
		{
			resId = pStem->getResidue(meConnection);
		}
		return resId;
	}

	bool StemConnection::isValid() const
	{
		return (NULL != mpStructure && meConnection != Stem::eUNDEFINED_CONNECTION);
	}

	bool StemConnection::operator== (const StemConnection &other) const
	{
		bool bEqual = (mpStructure == other.mpStructure)
			&& (meConnection == other.meConnection);
		return bEqual;
	}
};

